function mpc = delta_theta_bounds_0
mpc.version = '2';
mpc.baseMVA = 100.0;

mpc.bus = [
	1    3    0    0    0    0    1    1    0    240    1    1.1    0.9;
	2    1    100  100  0    0    1    1    0    240    1    1.1    0.9;
];

mpc.gen = [
	1 100  100  300 -300 1    100    1    300  -300;
];

mpc.gencost = [
	   2    0    0    3    1    1    1;
];

mpc.branch = [
	   1    2    0.1    0.1    0.1    500  500  500  0    0    1    0    0;
];
