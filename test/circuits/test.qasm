
OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];

cx q[0],q[1];
cx q[1],q[2];
cx q[2],q[1];
cx q[1],q[0];
u1(-0.7853981633974483) q[0];
cx q[1],q[0];
u2(0.7853981633974483,3.141592653589793) q[0];
cx q[1],q[0];
cx q[0],q[1];
u1(7.0685834705770345) q[2];
