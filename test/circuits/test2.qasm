// i 4 3 2 1 0
// o 1 4 3 0 2
OPENQASM 2.0;
include "qelib1.inc";
qreg q[5];
creg c[5];
cx q[4],q[3];
cx q[3],q[4];
cx q[4],q[3];
u3(3.141592653589793,0.0,3.141592653589793) q[4];
cx q[3],q[4];
cx q[4],q[3];
cx q[3],q[4];
cx q[1],q[3];
u3(3.141592653589793,0.0,3.141592653589793) q[1];
cx q[1],q[2];
cx q[2],q[1];
cx q[1],q[2];
u2(0.0,3.141592653589793) q[3];
cx q[4],q[3];
u1(-0.7853981633974483) q[3];
cx q[1],q[3];
u1(0.7853981633974483) q[3];
cx q[4],q[3];
u1(-0.7853981633974483) q[3];
cx q[3],q[1];
cx q[1],q[3];
u2(0.0,3.9269908169872414) q[1];
cx q[1],q[2];
u1(0.7853981633974483) q[4];
cx q[3],q[4];
u1(0.7853981633974483) q[3];
u1(-0.7853981633974483) q[4];
cx q[3],q[4];