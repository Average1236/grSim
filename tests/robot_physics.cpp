#include "robot.h"
#include <QApplication>
#include <QTemporaryDir>
#include <cmath>
#include <iostream>
#include <stdexcept>

static void near(double actual, double expected, double tolerance, const char* label)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::cerr << label << ": " << actual << " != " << expected << '\n';
        throw std::runtime_error(label);
    }
}

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QTemporaryDir runtime;
    if (argc != 2 || !runtime.isValid()) return 1;
    qputenv("KUN2_RUNTIME_DIR", runtime.path().toUtf8());
    qputenv("KUN2_GRSIM_CONFIG_DIR", argv[1]);
    ConfigWidget cfg;
    cfg.loadRobotSettings("Parsian");
    PWorld world(0.001, 0, nullptr, 2);
    PBall ball(2, 0, 0.1, 0.0215, 0.043, 1, 1, 1);
    world.addObject(&ball);
    Robot robot(&world, &ball, &cfg, 0, 0, ROBOT_START_Z((&cfg)), 1, 1, 1, 1, -1, 1);

    // Inspect actual ODE bodies, not just values copied into the configuration.
    const dBodyID bodies[] = {robot.chassis->body, robot.dummy->body, robot.kicker->box->body,
        robot.wheels[0]->cyl->body, robot.wheels[1]->cyl->body,
        robot.wheels[2]->cyl->body, robot.wheels[3]->cyl->body};
    double massSum = 0, mx = 0, my = 0, yawAboutOrigin = 0;
    for (dBodyID body : bodies) {
        dMass mass;
        dBodyGetMass(body, &mass);
        const dReal* pos = dBodyGetPosition(body);
        const dReal* rotation = dBodyGetRotation(body);
        // Independent evaluation of e_z^T R I R^T e_z + parallel-axis term.
        double yaw = 0;
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                yaw += rotation[8+i] * mass.I[4*i+j] * rotation[8+j];
        massSum += mass.mass;
        mx += mass.mass * pos[0];
        my += mass.mass * pos[1];
        yawAboutOrigin += yaw + mass.mass * (pos[0]*pos[0] + pos[1]*pos[1]);
    }
    near(massSum, 3.3, 1e-7, "assembled mass kg");
    near(yawAboutOrigin - (mx*mx + my*my)/massSum, 0.008, 1e-10, "assembled yaw inertia");
    for (auto* wheel : robot.wheels) {
        const dReal* pos = dBodyGetPosition(wheel->cyl->body);
        near(std::hypot(pos[0], pos[1]), 0.0785, 1e-10, "wheel centre radius");
        dReal radius, length;
        dGeomCylinderGetParams(wheel->cyl->geom, &radius, &length);
        near(radius, 0.033, 1e-10, "wheel radius");
        dMass mass;
        dBodyGetMass(wheel->cyl->body, &mass);
        near(mass.I[10], 5e-5, 1e-12, "wheel axial inertia");
    }

    // A second profile must not change an already-created robot's motor/kinematics.
    cfg.loadRobotSettings("RoboIME2012");
    Robot other(&world, &ball, &cfg, 1, 0, ROBOT_START_Z((&cfg)), 1, 1, 1, 2, -1, 1);
    robot.setSpeed(0, 0, 1, false, 0);
    for (auto* wheel : robot.wheels) {
        wheel->step();
        near(dJointGetAMotorParam(wheel->motor, dParamVel), 0.0785/0.033, 1e-9, "yaw wheel speed");
        near(dJointGetAMotorParam(wheel->motor, dParamFMax), 0.3, 1e-10, "Parsian motor limit");
    }
    other.wheels[0]->step();
    near(dJointGetAMotorParam(other.wheels[0]->motor, dParamFMax), 0.2, 1e-10, "legacy motor limit");

    // Disable the servo and pin the chassis: an unpowered spinning wheel must
    // dissipate energy according to I*d(omega)/dt = -b*omega, in both directions.
    dBodySetKinematic(robot.chassis->body);
    world.initAllObjects();
    auto* wheel = robot.wheels[0];
    for (auto* w : robot.wheels) dJointSetAMotorParam(w->motor, dParamFMax, 0);
    dVector3 axis;
    dJointGetHingeAxis(wheel->joint, axis);
    for (double initial : {-40.0, 40.0}) {
        dBodySetAngularVel(wheel->cyl->body, initial*axis[0], initial*axis[1], initial*axis[2]);
        const double before = dJointGetHingeAngleRate(wheel->joint);
        for (int step=0; step<100; ++step) {
            wheel->applyDamping();
            world.step(0.001);
        }
        const double after = dJointGetHingeAngleRate(wheel->joint);
        near(after/before, std::exp(-6.5e-5/5e-5*0.1), 0.003, "unpowered wheel decay");
    }
    std::cout << "PASS: ODE mass, yaw/wheel inertia, geometry, profile isolation, wheel decay\n";
}
