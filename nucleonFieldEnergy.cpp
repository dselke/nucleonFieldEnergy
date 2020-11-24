/* nucleonFieldEnergy.cpp : win32 console application written in Visual Studio 2019. *//* David L. Selke 2014 */
/* the program ELLIPTIC.C by E. Dennison has been included with original comments below.*/
/* it is used for magnetic field computations, specifically the complete elliptic integrals of 1st and 2nd kind */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define SPEED_OF_LIGHT 299792458.0
#define ELECTRON_CHARGE 1.602176565E-19
#define PLANCKS_CONSTANT 6.62606957e-34
#define EPSILON_NOUGHT 8.85418781e-12
#define COULOMBS_CONSTANT_K 8.987552e9
#define MU_0 1.25663706212e-6
#define PI 3.141592
#define NEUTRON_RADIUS (8.775e-16)*(1.410606743/1.405220)
#define PROTON_RADIUS (8.775e-16)*(1.410606743/1.405220)
//README:  following three #defines are the resolution knobs of the sim.
//PROTON_RADII_IN_SIM defines how far apart the edges of the simulated area are.  The orbitspheres are centered in the middle.
//num_points_grid_line says how many grid cubes tall the proton is (diameter)
//num_points_circle says how many points of charge are equally spaced around each circle on the orbitsphere surface
#define PROTON_RADII_IN_SIM 5
#define num_points_grid_line 16
#define num_points_circle 16
//
int num_points_sphere = 2 * (num_points_circle / 2 - 1) * (num_points_circle / 2);
double grid_step = 2 * PROTON_RADIUS / (double)num_points_grid_line;

typedef struct vector
{
	double magnitude_x;
	double magnitude_y;
	double magnitude_z;
	double magnitude;
	double direction_theta;
	double direction_phi;
} vector;

typedef struct point
{
	double x;
	double y;
	double z;
} point;

typedef struct orbitsphere
{
	double center_x;
	double center_y;
	double center_z;
	double radius;
	double rotation_theta;
	double rotation_phi;
	double velocity_angle_from_radius;
	double* points_x;
	double* points_y;
	double* points_z;
	double* charges;
	double* masses;
	double* angle_theta;
	double* angle_phi;
} orbitsphere;

void computePointsNeutron(orbitsphere* orbit);
void computePointsNeutronQuark(orbitsphere* orbit);
void computePointsProton(orbitsphere* orbit);
void computePointsProtonQuark(orbitsphere* orbit);
void computeChargeWithAreaUniform(orbitsphere* orbit);
void computeChargeWithAreaProton(orbitsphere* orbit);
void computeChargeWithAreaNeutron(orbitsphere* orbit);
void computeChargeWithAreaUpQuark(orbitsphere* orbit);
void computeChargeWithAreaUpQuark2(orbitsphere* orbit);
void computeChargeWithAreaDownQuark(orbitsphere* orbit);
void computeChargeWithAreaDownQuark2(orbitsphere* orbit);
double computeBFieldEnergyAtPointProton(orbitsphere* a, double x, double y, double z);
void computeBFieldVectorAtPointProton(orbitsphere* o1, double x1, double y1, double z1, vector* result);
double computeEFieldEnergyAtPointProton(orbitsphere* a, double x, double y, double z);
void computeEFieldVectorAtPointProton(orbitsphere* a, double x, double y, double z, vector* result);
double computeBFieldEnergyAtPointNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c, double x, double y, double z);
void computeBFieldVectorAtPointNeutron(orbitsphere* o1, orbitsphere* o2, orbitsphere* o3, double x1, double y1, double z1, vector* result);
double computeEFieldEnergyAtPointNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c, double x, double y, double z);
void computeEFieldVectorAtPointNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c, double x, double y, double z, vector* result);
double computeBFieldEnergyProton(orbitsphere* a);
double computeBFieldEnergyNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c);
double computeEFieldEnergyProton(orbitsphere* a);
double computeEFieldEnergyNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c);
void setCenter(orbitsphere* orbit, double x, double y, double z);
void setRotation(orbitsphere* orbit, double theta, double phi);
void setRotationDegrees(orbitsphere* orbit, double theta, double phi);
void crossProduct(vector* a, vector* b, vector* result);
void rotate2D(vector* a, double angle, vector* result);
void rotate3D(vector* vec, double angle, vector* axis_unit, point* axis_point, vector* result);
void setVectorMagnitude(vector* a, double magnitude);
void setVectorDirection(vector* a, double theta, double phi);
void vectorAdd(vector* a, vector* b, vector* result);
double completeEllipticIntegralFirstKind(double k);
double completeEllipticIntegralSecondKind(double k);
double drf(double x, double y, double z, int* piErr);
double drd(double x, double y, double z, int* piErr);

int main(void)
{
	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double degrees = 0;
	double radians = 0;
	double energy;
	vector* net_torque = new vector;
	//vector* moment_component;
	orbitsphere* a = new orbitsphere;
	a->points_x = (double*)malloc(num_points_sphere * sizeof(double));
	a->points_y = (double*)malloc(num_points_sphere * sizeof(double));
	a->points_z = (double*)malloc(num_points_sphere * sizeof(double));
	a->charges = (double*)malloc(num_points_sphere * sizeof(double));
	a->angle_theta = (double*)malloc(num_points_sphere * sizeof(double));
	a->angle_phi = (double*)malloc(num_points_sphere * sizeof(double));
	setRotation(a, 0, 0);
	setCenter(a, 0, 0, 0);
	printf("num_points_sphere %d\n", num_points_sphere);
	printf("using proton radius %e\n", PROTON_RADIUS);
	printf("using neutron radius %e\n", NEUTRON_RADIUS);

	printf("PROTON_RADII_IN_SIM %d\n", PROTON_RADII_IN_SIM);
	printf("num_points_grid_line %d\n", num_points_grid_line);

	int num_points_sphere = 2 * (num_points_circle / 2 - 1) * (num_points_circle / 2);
	int num_points_grid = num_points_grid_line * num_points_grid_line * num_points_grid_line;
	double grid_step = 2 * PROTON_RADIUS / (double)num_points_grid_line;


	computePointsProton(a);
	computeChargeWithAreaProton(a);
	a->velocity_angle_from_radius = PI / 2;
	energy = computeEFieldEnergyProton(a);
	printf("Proton E field energy within 95 percent of radius from center or outside 105 percent\n %e (MKS units)\n", energy);
	computeBFieldEnergyProton(a);

	energy = computeBFieldEnergyProton(a);
	printf("Proton B field energy at points at least 1/32 grid step from all current loop planes\n %e (MKS units)\n", energy);
	computeBFieldEnergyProton(a);
	orbitsphere* u = new orbitsphere;//this is not needed; a uniform distribution of charge on the sphere for test purposes
	//u->points_x = (double*)malloc(num_points_sphere * sizeof(double));
	//u->points_y = (double*)malloc(num_points_sphere * sizeof(double));
	//u->points_z = (double*)malloc(num_points_sphere * sizeof(double));
	//u->charges = (double*)malloc(num_points_sphere * sizeof(double));
	//u->angle_theta = (double*)malloc(num_points_sphere * sizeof(double));
	//u->angle_phi = (double*)malloc(num_points_sphere * sizeof(double));
	//computePointsProton(u);
	//computeChargeWithAreaUniform(u);
	//energy = computeEFieldEnergyProton(u);
	//printf("Uniform sphere E field energy within 95 percent of radius from center or outside 105 percent\n %e (MKS units)\n", energy);

	orbitsphere* up = new orbitsphere;
	up->points_x = (double*)malloc(num_points_sphere * sizeof(double));
	up->points_y = (double*)malloc(num_points_sphere * sizeof(double));
	up->points_z = (double*)malloc(num_points_sphere * sizeof(double));
	up->charges = (double*)malloc(num_points_sphere * sizeof(double));
	up->angle_theta = (double*)malloc(num_points_sphere * sizeof(double));
	up->angle_phi = (double*)malloc(num_points_sphere * sizeof(double));
	setRotation(up, 0, 0);
	setCenter(up, 0, 0, 0);
	computePointsNeutronQuark(up);
	computeChargeWithAreaUpQuark(up);
	up->velocity_angle_from_radius = PI / 2;
	energy = 0;
	orbitsphere* down1 = new orbitsphere;
	down1->points_x = (double*)malloc(num_points_sphere * sizeof(double));
	down1->points_y = (double*)malloc(num_points_sphere * sizeof(double));
	down1->points_z = (double*)malloc(num_points_sphere * sizeof(double));
	down1->charges = (double*)malloc(num_points_sphere * sizeof(double));
	down1->angle_theta = (double*)malloc(num_points_sphere * sizeof(double));
	down1->angle_phi = (double*)malloc(num_points_sphere * sizeof(double));
	setRotation(down1, 0, 0);
	setCenter(down1, 0, 0, 0);
	computePointsNeutronQuark(down1);
	computeChargeWithAreaDownQuark(down1);
	down1->velocity_angle_from_radius = PI / 2;
	orbitsphere* down2 = new orbitsphere;
	down2->points_x = (double*)malloc(num_points_sphere * sizeof(double));
	down2->points_y = (double*)malloc(num_points_sphere * sizeof(double));
	down2->points_z = (double*)malloc(num_points_sphere * sizeof(double));
	down2->charges = (double*)malloc(num_points_sphere * sizeof(double));
	down2->angle_theta = (double*)malloc(num_points_sphere * sizeof(double));
	down2->angle_phi = (double*)malloc(num_points_sphere * sizeof(double));
	setRotation(down2, 0, 0);
	setCenter(down2, 0, 0, 0);
	computePointsNeutronQuark(down2);
	computeChargeWithAreaDownQuark2(down2);
	down2->velocity_angle_from_radius = 3 * PI / 2;
	energy = computeEFieldEnergyNeutron(up, down1, down2);
	printf("Neutron E field energy within 95 percent of radius from center or outside 105 percent\n %e (MKS units)\n", energy);
	energy = computeBFieldEnergyNeutron(up, down1, down2);
	printf("Neutron B field energy at points at least 1/32 grid step from all current loop planes\n %e (MKS units)\n", energy);//[ ]TODO make quark rotation be respected in B field

	
	return 0;
}

void computePointsNeutron(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double circle_rad_step = (double)2 * PI / num_points_circle;
	double theta;
	double phi;
	orbit->radius = NEUTRON_RADIUS;
	theta = 0;
	phi = 0;
	for (i = 1; i < num_points_circle / 2; i++)
	{
		for (j = 0; j < num_points_circle; j++)
		{
			theta = i * circle_rad_step;
			phi = j * circle_rad_step;
			orbit->points_x[i - 1 + j * ((num_points_circle / 2) - 1)] = orbit->center_x + orbit->radius * sin(theta + orbit->rotation_theta) * cos(phi + orbit->rotation_phi);
			orbit->points_y[i - 1 + j * ((num_points_circle / 2) - 1)] = orbit->center_y + orbit->radius * sin(theta + orbit->rotation_theta) * sin(phi + orbit->rotation_phi);
			orbit->points_z[i - 1 + j * ((num_points_circle / 2) - 1)] = orbit->center_z + orbit->radius * cos(theta + orbit->rotation_theta);
			orbit->angle_theta[i - 1 + j * ((num_points_circle / 2) - 1)] = theta;
			orbit->angle_phi[i - 1 + j * ((num_points_circle / 2) - 1)] = phi;
		}
	}

	return;
}

void computePointsProton(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	int k = 0;
	int points_halfway = 0;
	int last_point = 0;
	double circle_rad_step = (double)2 * PI / num_points_circle;
	double theta;
	double phi;
	orbit->radius = PROTON_RADIUS;
	theta = 0;
	phi = 0;
	for (i = 1; i < num_points_circle / 2; i++)
	{
		for (j = 0; j < num_points_circle; j++)
		{
			theta = i * circle_rad_step;
			phi = j * circle_rad_step;
			orbit->points_x[i - 1 + j * ((num_points_circle / 2) - 1)] = orbit->center_x + orbit->radius * sin(theta + orbit->rotation_theta) * cos(phi + orbit->rotation_phi);
			orbit->points_y[i - 1 + j * ((num_points_circle / 2) - 1)] = orbit->center_y + orbit->radius * sin(theta + orbit->rotation_theta) * sin(phi + orbit->rotation_phi);
			orbit->points_z[i - 1 + j * ((num_points_circle / 2) - 1)] = orbit->center_z + orbit->radius * cos(theta + orbit->rotation_theta);
			orbit->angle_theta[i - 1 + j * ((num_points_circle / 2) - 1)] = theta;
			orbit->angle_phi[i - 1 + j * ((num_points_circle / 2) - 1)] = phi;
		}
	}
	return;
}

void computeChargeWithAreaNeutron(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;

	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double temp;
	double temp2 = 0;
	double totalCharge = 0;
	double totalChargePlus = 0;
	double totalChargeMinus = 0;
	double relativeChargePlus = 0;
	double relativeChargeMinus = 0;
	double* chargesMinus = (double*)malloc(num_points_sphere * sizeof(double));
	double area;
	double zone_h1;
	double zone_h2;
	double zone_height;
	double particle_radius = sqrt(orbit->points_x[0] * orbit->points_x[0] + orbit->points_y[0] * orbit->points_y[0] + orbit->points_z[0] * orbit->points_z[0]);
	for (i = 0; i < num_points_sphere; i++)
	{
		zone_h1 = particle_radius * cos(orbit->angle_theta[i] - circle_rad_step / 2);
		zone_h2 = particle_radius * cos(orbit->angle_theta[i] + circle_rad_step / 2);
		/*z coordinates of the band that reaches halfway to the next higher and to the next lower point in elevation *//*for computing the area represented by each point */
		if (zone_h1 > zone_h2)
		{
			zone_height = zone_h1 - zone_h2;
		}
		else
		{
			zone_height = zone_h2 - zone_h1;
		}
		area = 2 * PI * particle_radius * zone_height / num_points_circle;
		temp = (2.0 / 3.0) * (1 + sin(orbit->angle_theta[i] + orbit->rotation_theta) * sin(orbit->angle_phi[i] + orbit->rotation_phi));
		temp = temp * area;
		relativeChargePlus += temp;
		orbit->charges[i] = temp;
		temp2 = (-1.0 / 3.0) * (1 + sin(orbit->angle_theta[i] + orbit->rotation_theta) * cos(orbit->angle_phi[i] + orbit->rotation_phi));
		temp2 -= (1.0 / 3.0) * (1 + cos(orbit->angle_theta[i] + orbit->rotation_theta));
		temp2 = temp2 * area;
		relativeChargeMinus += temp2;
		chargesMinus[i] = temp2;
	}
	i = 0;
	for (i = 0; i < num_points_sphere; i++)
	{
		orbit->charges[i] = (orbit->charges[i] * (2.0 / 3.0) * ELECTRON_CHARGE / relativeChargePlus);
		totalChargePlus += orbit->charges[i];
		totalChargeMinus -= (chargesMinus[i] * (2.0 / 3.0) * ELECTRON_CHARGE / relativeChargeMinus);
		orbit->charges[i] -= (chargesMinus[i] * (2.0 / 3.0) * ELECTRON_CHARGE / relativeChargeMinus);
		totalCharge += orbit->charges[i];
	}
	free(chargesMinus);
	return;
}

void computeChargeWithAreaProton(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double temp;
	double totalCharge = 0;
	double relativeCharge = 0;
	double area;
	/* the space between points is less at the poles: multiply each point by the "square" area that it is at the center of */double zone_h1;
	double zone_h2;
	double zone_height;
	double particle_radius = sqrt(orbit->points_x[0] * orbit->points_x[0] + orbit->points_y[0] * orbit->points_y[0] + orbit->points_z[0] * orbit->points_z[0]);
	for (i = 0; i < num_points_sphere; i++)
	{
		zone_h1 = particle_radius * cos(orbit->angle_theta[i] - circle_rad_step / 2);
		zone_h2 = particle_radius * cos(orbit->angle_theta[i] + circle_rad_step / 2);
		/*z coordinates of the band that reaches halfway to the next higher and to the next lower point in elevation *//*for computing the area represented by each point */
		if (zone_h1 > zone_h2)
		{
			zone_height = zone_h1 - zone_h2;
		}
		else
		{
			zone_height = zone_h2 - zone_h1;
		}
		area = 2 * PI * particle_radius * zone_height / num_points_circle;
		temp = (2.0 / 3.0) * (1 + sin(orbit->angle_theta[i] + orbit->rotation_theta) * sin(orbit->angle_phi[i] + orbit->rotation_phi));
		temp += (2.0 / 3.0) * (1 + sin(orbit->angle_theta[i] + orbit->rotation_theta) * cos(orbit->angle_phi[i] + orbit->rotation_phi));
		temp -= (1.0 / 3.0) * (1 + cos(orbit->angle_theta[i] + orbit->rotation_theta));
		temp = temp * area;
		orbit->charges[i] = temp;
		relativeCharge += temp;
	}
	i = 0;
	for (i = 0; i < num_points_sphere; i++)
	{
		orbit->charges[i] = orbit->charges[i] * ELECTRON_CHARGE / relativeCharge;
		totalCharge += orbit->charges[i];
	}
	return;
}

void computeChargeWithAreaUniform(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double temp;
	double totalCharge = 0;
	double relativeCharge = 0;
	double area;
	/* the space between points is less at the poles: multiply each point by the "square" area that it is at the center of */double zone_h1;
	double zone_h2;
	double zone_height;
	double particle_radius = sqrt(orbit->points_x[0] * orbit->points_x[0] + orbit->points_y[0] * orbit->points_y[0] + orbit->points_z[0] * orbit->points_z[0]);
	for (i = 0; i < num_points_sphere; i++)
	{
		zone_h1 = particle_radius * cos(orbit->angle_theta[i] - circle_rad_step / 2);
		zone_h2 = particle_radius * cos(orbit->angle_theta[i] + circle_rad_step / 2);
		/*z coordinates of the band that reaches halfway to the next higher and to the next lower point in elevation *//*for computing the area represented by each point */
		if (zone_h1 > zone_h2)
		{
			zone_height = zone_h1 - zone_h2;
		}
		else
		{
			zone_height = zone_h2 - zone_h1;
		}
		area = 2 * PI * particle_radius * zone_height / num_points_circle;
		temp = 1;
		temp = temp * area;
		orbit->charges[i] = temp;
		relativeCharge += temp;
	}
	i = 0;
	for (i = 0; i < num_points_sphere; i++)
	{
		orbit->charges[i] = orbit->charges[i] * ELECTRON_CHARGE / relativeCharge;
		totalCharge += orbit->charges[i];
	}
	printf("uniform total charge was: %e\n", totalCharge);
	return;
}


void computeChargeWithAreaUpQuark(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double temp;
	double totalCharge = 0;
	double relativeCharge = 0;
	double area;
	double zone_h1;
	double zone_h2;
	double zone_height;
	double particle_radius = sqrt(orbit->points_x[0] * orbit->points_x[0] + orbit->points_y[0] * orbit->points_y[0] + orbit->points_z[0] * orbit->points_z[0]);

	for (i = 0; i < num_points_sphere; i++)
	{
		zone_h1 = particle_radius * cos(orbit->angle_theta[i] - circle_rad_step / 2);
		zone_h2 = particle_radius * cos(orbit->angle_theta[i] + circle_rad_step / 2);
		if (zone_h1 > zone_h2)
		{
			zone_height = zone_h1 - zone_h2;
		}
		else
		{
			zone_height = zone_h2 - zone_h1;
		}
		area = 2 * PI * particle_radius * zone_height / num_points_circle;
		temp = (2.0 / 3.0) * (1 + sin(orbit->angle_theta[i] + orbit->rotation_theta) * sin(orbit->angle_phi[i] + orbit->rotation_phi));
		temp = temp * area;
		orbit->charges[i] = temp;
		relativeCharge += temp;
	}
	i = 0;
	for (i = 0; i < num_points_sphere; i++)
	{
		orbit->charges[i] = orbit->charges[i] * (2.0 / 3.0) * ELECTRON_CHARGE / relativeCharge;
		totalCharge += orbit->charges[i];
	}
	return;
}

void computeChargeWithAreaUpQuark2(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double temp;
	double totalCharge = 0;
	double relativeCharge = 0;
	double area;
	double zone_h1;
	double zone_h2;
	double zone_height;
	double particle_radius = sqrt(orbit->points_x[0] * orbit->points_x[0] + orbit->points_y[0] * orbit->points_y[0] + orbit->points_z[0] * orbit->points_z[0]);
	for (i = 0; i < num_points_sphere; i++)
	{
		zone_h1 = particle_radius * cos(orbit->angle_theta[i] - circle_rad_step / 2);
		zone_h2 = particle_radius * cos(orbit->angle_theta[i] + circle_rad_step / 2);
		if (zone_h1 > zone_h2) {
			zone_height = zone_h1 - zone_h2;
		}
		else
		{
			zone_height = zone_h2 - zone_h1;
		}
		area = 2 * PI * particle_radius * zone_height / num_points_circle;
		temp = (2.0 / 3.0) * (1 + sin(orbit->angle_theta[i] + orbit->rotation_theta) * cos(orbit->angle_phi[i] + orbit->rotation_phi));
		temp = temp * area;
		orbit->charges[i] = temp;
		relativeCharge += temp;
	}
	i = 0;
	for (i = 0; i < num_points_sphere; i++)
	{
		orbit->charges[i] = orbit->charges[i] * (2.0 / 3.0) * ELECTRON_CHARGE / relativeCharge;
		totalCharge += orbit->charges[i];
	}
	return;
}

void computeChargeWithAreaDownQuark(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double temp;
	double totalCharge = 0;
	double relativeCharge = 0;
	double area;
	double zone_h1;
	double zone_h2;
	double zone_height;
	double particle_radius = sqrt(orbit->points_x[0] * orbit->points_x[0] + orbit->points_y[0] * orbit->points_y[0] + orbit->points_z[0] * orbit->points_z[0]);
	for (i = 0; i < num_points_sphere; i++)
	{
		zone_h1 = particle_radius * cos(orbit->angle_theta[i] - circle_rad_step / 2);
		zone_h2 = particle_radius * cos(orbit->angle_theta[i] + circle_rad_step / 2);
		if (zone_h1 > zone_h2)
		{
			zone_height = zone_h1 - zone_h2;
		}
		else
		{
			zone_height = zone_h2 - zone_h1;
		}
		area = 2 * PI * particle_radius * zone_height / num_points_circle;
		temp = (-1.0 / 3.0) * (1 + sin(orbit->angle_theta[i] + orbit->rotation_theta) * cos(orbit->angle_phi[i] + orbit->rotation_phi));
		temp = temp * area;
		orbit->charges[i] = temp;
		relativeCharge += temp;
	}
	i = 0;
	for (i = 0; i < num_points_sphere; i++)
	{
		orbit->charges[i] = orbit->charges[i] * (1.0 / 3.0) * ELECTRON_CHARGE / relativeCharge;
		totalCharge += orbit->charges[i];
	}
	return;
}

void computeChargeWithAreaDownQuark2(orbitsphere* orbit)
{
	int i = 0;
	int j = 0;
	double circle_rad_step = (double)2 * 3.141592 / num_points_circle;
	double temp;
	double totalCharge = 0;
	double relativeCharge = 0;
	double area;
	double zone_h1;
	double zone_h2;
	double zone_height;
	double particle_radius = sqrt(orbit->points_x[0] * orbit->points_x[0] + orbit->points_y[0] * orbit->points_y[0] + orbit->points_z[0] * orbit->points_z[0]);
	for (i = 0; i < num_points_sphere; i++)
	{
		zone_h1 = particle_radius * cos(orbit->angle_theta[i] - circle_rad_step / 2);
		zone_h2 = particle_radius * cos(orbit->angle_theta[i] + circle_rad_step / 2);
		if (zone_h1 > zone_h2)
		{
			zone_height = zone_h1 - zone_h2;
		}
		else
		{
			zone_height = zone_h2 - zone_h1;
		}
		area = 2 * PI * particle_radius * zone_height / num_points_circle;
		temp = (-1.0 / 3.0) * (1 + cos(orbit->angle_theta[i] + orbit->rotation_theta));
		temp = temp * area;
		orbit->charges[i] = temp;
		relativeCharge += temp;
	}
	i = 0;
	for (i = 0; i < num_points_sphere; i++)
	{
		orbit->charges[i] = orbit->charges[i] * (1.0 / 3.0) * ELECTRON_CHARGE / relativeCharge;
		totalCharge += orbit->charges[i];
	}
	return;
}

void setCenter(orbitsphere* orbit, double x, double y, double z)
{
	orbit->center_x = x;
	orbit->center_y = y;
	orbit->center_z = z;
	return;
}

void setRotation(orbitsphere* orbit, double theta, double phi) {
	orbit->rotation_theta = theta;
	orbit->rotation_phi = phi;
	return;
}

void setRotationDegrees(orbitsphere* orbit, double theta, double phi)
{
	orbit->rotation_theta = (360 * theta) / (2 * PI);
	orbit->rotation_phi = (360 * phi) / (2 * PI);
	return;
}

void crossProduct(vector* a, vector* b, vector* result) {
	result->magnitude_x = a->magnitude_y * b->magnitude_z - a->magnitude_z * b->magnitude_y;
	result->magnitude_y = a->magnitude_z * b->magnitude_x - a->magnitude_x * b->magnitude_z;
	result->magnitude_z = a->magnitude_x * b->magnitude_y - a->magnitude_y * b->magnitude_x;
	return;
}

void rotate2D(vector* a, double angle, vector* result)
{
	double rotMatrix11 = cos(angle);
	double rotMatrix12 = -1.0 * sin(angle);
	double rotMatrix21 = sin(angle);
	double rotMatrix22 = cos(angle);
	result->magnitude_x = a->magnitude_x * rotMatrix11 + a->magnitude_y * rotMatrix12;
	result->magnitude_y = a->magnitude_x * rotMatrix21 + a->magnitude_y * rotMatrix22;
	result->magnitude_z = 0.0;
}

//formula from http://www.cs.berkeley.edu/~ug/slide/pipeline/assignments/as5/rotation.html
//p' = (p + (a x (a x p))) + cos(b) (-(a x (a x p))) + sin(b) (a x p) 
void rotate3D(vector* vec, double angle, vector* axis_unit, point* axis_point, vector* result)
{
	vector a_cross_p;
	vector a_cross_a_cross_p;
	double cos_b;
	double sin_b;
	vector temp;
	vector temp2;
	vector temp3;
	//TODO remove these
	//double dot_product;
	//double magnitude_product;
	//double check_angle;
	//
	//printf("in rotate3D angle is %e\n", angle);//TODO remove
	cos_b = cos(angle);
	sin_b = sin(angle);
	//first term
	crossProduct(axis_unit, vec, &a_cross_p);
	crossProduct(axis_unit, &a_cross_p, &a_cross_a_cross_p);
	vectorAdd(vec, &a_cross_a_cross_p, &temp);
	//second term
	temp2.magnitude_x = -a_cross_a_cross_p.magnitude_x * cos_b;
	temp2.magnitude_y = -a_cross_a_cross_p.magnitude_y * cos_b;
	temp2.magnitude_z = -a_cross_a_cross_p.magnitude_z * cos_b;
	//third term
	temp3.magnitude_x = a_cross_p.magnitude_x * sin_b;
	temp3.magnitude_y = a_cross_p.magnitude_y * sin_b;
	temp3.magnitude_z = a_cross_p.magnitude_z * sin_b;

	result->magnitude_x = temp.magnitude_x + temp2.magnitude_x + temp3.magnitude_x;
	result->magnitude_y = temp.magnitude_y + temp2.magnitude_y + temp3.magnitude_y;
	result->magnitude_z = temp.magnitude_z + temp2.magnitude_z + temp3.magnitude_z;
	//TODO remove
	//dot_product = dotProduct(vec, result);
	//magnitude_product = sqrt(vec->x*vec->x + vec->y*vec->y + vec->z*vec->z)*sqrt(result->x*result->x + result->y*result->y + result->z*result->z);
	//check_angle = acos(dot_product / magnitude_product);
	//printf("    check_angle is %e\n", check_angle);
	//
	return;
}

void setVectorMagnitude(vector* a, double magnitude)
{
	double starting_magnitude;
	double ratio;
	starting_magnitude = sqrt(a->magnitude_x * a->magnitude_x + a->magnitude_y * a->magnitude_y + a->magnitude_z * a->magnitude_z);
	ratio = magnitude / starting_magnitude;
	a->magnitude_x = a->magnitude_x * ratio;
	a->magnitude_y = a->magnitude_y * ratio;
	a->magnitude_z = a->magnitude_z * ratio;
	a->magnitude = magnitude;
}

void setVectorDirection(vector* a, double theta, double phi) {
	double magnitude;
	magnitude = sqrt(a->magnitude_x * a->magnitude_x + a->magnitude_y * a->magnitude_y + a->magnitude_z * a->magnitude_z);
	a->magnitude_z = magnitude * cos(phi);
	a->magnitude_x = magnitude * sin(phi) * cos(theta);
	a->magnitude_y = magnitude * sin(phi) * sin(theta);
}

void vectorAdd(vector* a, vector* b, vector* result)
{
	result->magnitude_x = a->magnitude_x + b->magnitude_x;
	result->magnitude_y = a->magnitude_y + b->magnitude_y;
	result->magnitude_z = a->magnitude_z + b->magnitude_z;
}

void computePointsNeutronQuark(orbitsphere* orbit)
{
	computePointsNeutron(orbit);
}

void computePointsProtonQuark(orbitsphere* orbit)
{
	computePointsProton(orbit);
}

void computeEFieldVectorAtPointProton(orbitsphere* a, double x, double y, double z, vector* result)
{
	vector distanceVector;
	vector fieldAtPoint;
	vector resultantField;
	resultantField.magnitude_x = 0;
	resultantField.magnitude_y = 0;
	resultantField.magnitude_z = 0;
	result->magnitude_x = 0;
	result->magnitude_y = 0;
	result->magnitude_z = 0;
	result->magnitude = 0;
	if ((sqrt(x * x + y * y + z * z) > PROTON_RADIUS * 0.95) && (sqrt(x * x + y * y + z * z) < PROTON_RADIUS * 1.05))
	{//skip a point if it is too close to the charges
		//result vector is zero
		return;
	}
	for (int i = 0; i < num_points_sphere; i++)
	{
		distanceVector.magnitude_x = x - a->points_x[i];
		distanceVector.magnitude_y = y - a->points_y[i];
		distanceVector.magnitude_z = z - a->points_z[i];
		distanceVector.magnitude = sqrt(distanceVector.magnitude_x * distanceVector.magnitude_x + distanceVector.magnitude_y * distanceVector.magnitude_y + distanceVector.magnitude_z * distanceVector.magnitude_z);
		fieldAtPoint.magnitude_x = distanceVector.magnitude_x;//start with the same vector (because it is the right direction)
		fieldAtPoint.magnitude_y = distanceVector.magnitude_y;
		fieldAtPoint.magnitude_z = distanceVector.magnitude_x;
		setVectorMagnitude(&fieldAtPoint, (COULOMBS_CONSTANT_K * a->charges[i]) / (distanceVector.magnitude * distanceVector.magnitude));
		vectorAdd(&fieldAtPoint, &resultantField, &resultantField);
	}
	//now all the E field contributions of the points have been added up
	resultantField.magnitude = sqrt(resultantField.magnitude_x * resultantField.magnitude_x + resultantField.magnitude_y * resultantField.magnitude_y + resultantField.magnitude_z * resultantField.magnitude_z);
	result->magnitude_x = resultantField.magnitude_x;
	result->magnitude_y = resultantField.magnitude_y;
	result->magnitude_z = resultantField.magnitude_z;
	result->magnitude = sqrt(result->magnitude_x* result->magnitude_x + result->magnitude_y * result->magnitude_y + result->magnitude_z * result->magnitude_z);
	
	return;
}

double computeEFieldEnergyAtPointProton(orbitsphere* a, double x, double y, double z)
{
	vector EField;
	double energy;
	computeEFieldVectorAtPointProton(a, x, y, z, &EField);
	energy = 0.5*EPSILON_NOUGHT* EField.magnitude* EField.magnitude;//energy in a cubic meter
	return energy*(grid_step*grid_step*grid_step);//energy in a volume element
}

double computeEFieldEnergyAtPointNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c, double x, double y, double z)
{
	vector EField;
	double energy;
	computeEFieldVectorAtPointNeutron(a, b, c, x, y, z, &EField);
	energy = 0.5 * EPSILON_NOUGHT * EField.magnitude * EField.magnitude;//energy in a cubic meter
	return energy * (grid_step * grid_step * grid_step);//energy in a volume element
}


void computeEFieldVectorAtPointNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c, double x, double y, double z, vector* result)
{
	vector distanceVector;
	vector fieldAtPoint;
	vector resultantField;
	resultantField.magnitude_x = 0;
	resultantField.magnitude_y = 0;
	resultantField.magnitude_z = 0;
	result->magnitude_x = 0;
	result->magnitude_y = 0;
	result->magnitude_z = 0;
	result->magnitude = 0;
	if ((sqrt(x * x + y * y + z * z) > PROTON_RADIUS * 0.95) && (sqrt(x * x + y * y + z * z) < PROTON_RADIUS * 1.05))
	{//skip a point if it is too close to the charges
		//result vector is zero
		return;
	}
	for (int i = 0; i < num_points_sphere; i++)
	{
		distanceVector.magnitude_x = x - a->points_x[i];
		distanceVector.magnitude_y = y - a->points_y[i];
		distanceVector.magnitude_z = z - a->points_z[i];
		distanceVector.magnitude = sqrt(distanceVector.magnitude_x * distanceVector.magnitude_x + distanceVector.magnitude_y * distanceVector.magnitude_y + distanceVector.magnitude_z * distanceVector.magnitude_z);
		fieldAtPoint.magnitude_x = distanceVector.magnitude_x;//start with the same vector (because it is the right direction)
		fieldAtPoint.magnitude_y = distanceVector.magnitude_y;
		fieldAtPoint.magnitude_z = distanceVector.magnitude_x;
		setVectorMagnitude(&fieldAtPoint, (COULOMBS_CONSTANT_K * a->charges[i]) / (distanceVector.magnitude * distanceVector.magnitude));
		vectorAdd(&fieldAtPoint, &resultantField, &resultantField);

		distanceVector.magnitude_x = x - b->points_x[i];
		distanceVector.magnitude_y = y - b->points_y[i];
		distanceVector.magnitude_z = z - b->points_z[i];
		distanceVector.magnitude = sqrt(distanceVector.magnitude_x * distanceVector.magnitude_x + distanceVector.magnitude_y * distanceVector.magnitude_y + distanceVector.magnitude_z * distanceVector.magnitude_z);
		fieldAtPoint.magnitude_x = distanceVector.magnitude_x;//start with the same vector (because it is the right direction)
		fieldAtPoint.magnitude_y = distanceVector.magnitude_y;
		fieldAtPoint.magnitude_z = distanceVector.magnitude_x;
		setVectorMagnitude(&fieldAtPoint, (COULOMBS_CONSTANT_K * b->charges[i]) / (distanceVector.magnitude * distanceVector.magnitude));
		vectorAdd(&fieldAtPoint, &resultantField, &resultantField);

		distanceVector.magnitude_x = x - c->points_x[i];
		distanceVector.magnitude_y = y - c->points_y[i];
		distanceVector.magnitude_z = z - c->points_z[i];
		distanceVector.magnitude = sqrt(distanceVector.magnitude_x * distanceVector.magnitude_x + distanceVector.magnitude_y * distanceVector.magnitude_y + distanceVector.magnitude_z * distanceVector.magnitude_z);
		fieldAtPoint.magnitude_x = distanceVector.magnitude_x;//start with the same vector (because it is the right direction)
		fieldAtPoint.magnitude_y = distanceVector.magnitude_y;
		fieldAtPoint.magnitude_z = distanceVector.magnitude_x;
		setVectorMagnitude(&fieldAtPoint, (COULOMBS_CONSTANT_K * c->charges[i]) / (distanceVector.magnitude * distanceVector.magnitude));
		vectorAdd(&fieldAtPoint, &resultantField, &resultantField);

	}
	//now all the E field contributions of the points have been added up
	resultantField.magnitude = sqrt(resultantField.magnitude_x * resultantField.magnitude_x + resultantField.magnitude_y * resultantField.magnitude_y + resultantField.magnitude_z * resultantField.magnitude_z);
	result->magnitude_x = resultantField.magnitude_x;
	result->magnitude_y = resultantField.magnitude_y;
	result->magnitude_z = resultantField.magnitude_z;
	result->magnitude = sqrt(result->magnitude_x * result->magnitude_x + result->magnitude_y * result->magnitude_y + result->magnitude_z * result->magnitude_z);

	return;
}

double computeBFieldEnergyProton(orbitsphere* a)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double total_energy = 0;
	for (i = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i++)
	{
		for (j = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j++)
		{
			for (k = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k++)
			{
				total_energy += computeBFieldEnergyAtPointProton(a, (double)i * grid_step, (double)j * grid_step, (double)k * grid_step);
			}
		}
	}
	return total_energy;
}

double computeBFieldEnergyNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double total_energy = 0;
	for (i = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i++)
	{
		for (j = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j++)
		{
			for (k = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k++)
			{
				total_energy += computeBFieldEnergyAtPointNeutron(a, b, c, (double)i * grid_step, (double)j * grid_step, (double)k * grid_step);
			}
		}
	}
	return total_energy;
}


double computeEFieldEnergyProton(orbitsphere* a)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double total_energy = 0;
	for (i = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i++)
	{
		for (j = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j++)
		{
			for (k = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k++)
			{
				total_energy += computeEFieldEnergyAtPointProton(a, (double)i * grid_step, (double)j * grid_step, (double)k * grid_step);
			}
		}
	}
	return total_energy;
}

double computeEFieldEnergyNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c)
{
	int i = 0;
	int j = 0;
	int k = 0;
	double total_energy = 0;
	for (i = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); i++)
	{
		for (j = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); j++)
		{
			for (k = -PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k < PROTON_RADII_IN_SIM * (num_points_grid_line / 2); k++)
			{
				total_energy += computeEFieldEnergyAtPointNeutron(a, b, c, (double)i * grid_step, (double)j * grid_step, (double)k * grid_step);
			}
		}
	}
	return total_energy;
}

double computeBFieldEnergyAtPointProton(orbitsphere* a, double x, double y, double z)
{
	vector BField;
	double energy;
	computeBFieldVectorAtPointProton(a, x, y, z, &BField);
	energy = BField.magnitude * BField.magnitude/(2 * MU_0);//energy in a cubic meter
	return energy * (grid_step * grid_step * grid_step);//energy in a volume element
}
double computeBFieldEnergyAtPointNeutron(orbitsphere* a, orbitsphere* b, orbitsphere* c, double x, double y, double z)
{
	vector BField;
	double energy;
	computeBFieldVectorAtPointNeutron(a, b, c, x, y, z, &BField);
	energy = BField.magnitude * BField.magnitude / (2 * MU_0);//energy in a cubic meter
	return energy * (grid_step * grid_step * grid_step);//energy in a volume element
}
void computeBFieldVectorAtPointProton(orbitsphere* o1, double x1, double y1, double z1, vector* result)
{
	double fieldAtPointBxAxial;
	double fieldAtPointBrRadial;
	vector fieldAtPointXyz;
	vector resultantFieldXyz;
	vector loopCenter;
	vector loopCenterToFieldPoint;
	int pointIndexOne = 0;
	int pointIndexTwo = 0;
	double loopCharge = 0.0;
	double loopCurrent = 0.0;
	double revsPerSecond = 0;
	resultantFieldXyz.magnitude_x = 0;
	resultantFieldXyz.magnitude_y = 0;
	resultantFieldXyz.magnitude_z = 0;
	result->magnitude_x = 0;
	result->magnitude_y = 0;
	result->magnitude_z = 0;
	result->magnitude = 0;
	revsPerSecond = 1/(o1->radius*2*PI/SPEED_OF_LIGHT);//the angular velocity is same across orbitsphere

	for (int i = 0; i < num_points_sphere; i += num_points_circle)
	{
		pointIndexOne = i;
		pointIndexTwo = pointIndexOne + num_points_circle / 2;
		loopCenter.magnitude_x = (o1->points_x[pointIndexOne] + o1->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o1->points_y[pointIndexOne] + o1->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o1->points_z[pointIndexOne] + o1->points_z[pointIndexTwo]) / 2;
		loopCharge = 0.0;
		for (int j = 0; j < num_points_circle; j++)
		{
			loopCharge += o1->charges[i+j];
		}
		loopCurrent = loopCharge * revsPerSecond; //Coulombs per second (Amps)

	//The following code is from main() of ELLIPTIC.C.  There are modifications, but the unmodified code is at the end of this file in a comment.
	//	/* input parameters */
		loopCenter.magnitude_x = (o1->points_x[pointIndexOne] + o1->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o1->points_y[pointIndexOne] + o1->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o1->points_z[pointIndexOne] + o1->points_z[pointIndexTwo]) / 2;

		double temp = loopCenter.magnitude_x - o1->points_x[pointIndexOne];
		temp = temp * temp;
		double temp2 = loopCenter.magnitude_y - o1->points_y[pointIndexOne];
		temp2 = temp2 * temp2;
		double temp3 = loopCenter.magnitude_z - o1->points_z[pointIndexOne];
		temp3 = temp3 * temp3;
		temp = temp + temp2 + temp3;
		float a = (float)sqrt(temp);
		double axialOffset = o1->points_z[pointIndexOne];//z coord of the present loop[ ] TODO add axial offset back in?
		axialOffset = axialOffset - z1;//[ ] TODO check sign! field computed is relative to axialOffset
		float r = (float)sqrt(x1*x1 + y1*y1);
		if (axialOffset < grid_step/32)//if the plane of the field point is close to the plane of the loop: skip this point
		{//[ ] TODO tighten up this condition
			return;//result is 0 so far
		}
		loopCenterToFieldPoint.magnitude_x = x1 - loopCenter.magnitude_x;
		loopCenterToFieldPoint.magnitude_y = y1 - loopCenter.magnitude_y;
		double radius_angle = atan(loopCenterToFieldPoint.magnitude_x/ loopCenterToFieldPoint.magnitude_y);
		double alpha = r / a;
		double beta = axialOffset / a;
		double gamma = axialOffset / r;
		double al2 = alpha *alpha;
		double be2 = beta * beta;
		double q = pow(1 + alpha, 2) + be2;
		double rq = sqrt(q);
		double k = sqrt((alpha*4) / q);
		double fk = completeEllipticIntegralFirstKind(k);
		double ek = completeEllipticIntegralSecondKind(k);
		double Ht = loopCurrent / (2.0 * a * PI * rq);
		double Hx = Ht * (ek * (1 - al2 - be2) / (q - (alpha*4)) + fk);
		double Hr = (r == 0.0) ? (0.0) : (Ht * gamma * (ek * (1 + al2 + be2) / (q - (alpha * 4)) - fk));
		fieldAtPointBxAxial = Hx * MU_0;//convert from H to B field
		fieldAtPointBrRadial = Hr * MU_0;
		fieldAtPointXyz.magnitude_x = fieldAtPointBrRadial * sin(radius_angle);//[ ] TODO check on sin, cos
		fieldAtPointXyz.magnitude_y = fieldAtPointBrRadial * cos(radius_angle);//make sure not backwards
		fieldAtPointXyz.magnitude_z = fieldAtPointBxAxial;
		resultantFieldXyz.magnitude_x = resultantFieldXyz.magnitude_x + fieldAtPointXyz.magnitude_x;
		resultantFieldXyz.magnitude_y = resultantFieldXyz.magnitude_y + fieldAtPointXyz.magnitude_y;
		resultantFieldXyz.magnitude_z = resultantFieldXyz.magnitude_z + fieldAtPointXyz.magnitude_z;

	}
	result->magnitude_x = resultantFieldXyz.magnitude_x;
	result->magnitude_y = resultantFieldXyz.magnitude_y;
	result->magnitude_z = resultantFieldXyz.magnitude_z;
	result->magnitude = sqrt(result->magnitude_x* result->magnitude_x+ result->magnitude_y* result->magnitude_y+ result->magnitude_z* result->magnitude_z);
}
void computeBFieldVectorAtPointNeutron(orbitsphere* o1, orbitsphere* o2, orbitsphere* o3, double x1, double y1, double z1, vector* result) 
{
	double fieldAtPointBxAxial;
	double fieldAtPointBrRadial;
	vector fieldAtPointXyz;
	vector resultantFieldXyz;
	vector loopCenter;
	vector loopCenterToFieldPoint;
	int pointIndexOne = 0;
	int pointIndexTwo = 0;
	double loopCharge = 0.0;
	double loopCurrent = 0.0;
	double revsPerSecond = 0;
	resultantFieldXyz.magnitude_x = 0;
	resultantFieldXyz.magnitude_y = 0;
	resultantFieldXyz.magnitude_z = 0;
	result->magnitude_x = 0;
	result->magnitude_y = 0;
	result->magnitude_z = 0;
	result->magnitude = 0;
	revsPerSecond = 1 / (o1->radius * 2 * PI / SPEED_OF_LIGHT);//the angular velocity is same across orbitsphere
	//first quark's current loops
	for (int i = 0; i < num_points_sphere; i += num_points_circle)
	{
		pointIndexOne = i;
		pointIndexTwo = pointIndexOne + num_points_circle / 2;
		loopCenter.magnitude_x = (o1->points_x[pointIndexOne] + o1->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o1->points_y[pointIndexOne] + o1->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o1->points_z[pointIndexOne] + o1->points_z[pointIndexTwo]) / 2;
		loopCharge = 0.0;
		for (int j = 0; j < num_points_circle; j++)
		{
			loopCharge += o1->charges[i + j];
		}
		loopCurrent = loopCharge * revsPerSecond; //Coulombs per second (Amps)

	//The following code is from main() of ELLIPTIC.C.  There are modifications, but the unmodified code is at the end of this file in a comment.
	//	/* input parameters */
		loopCenter.magnitude_x = (o1->points_x[pointIndexOne] + o1->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o1->points_y[pointIndexOne] + o1->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o1->points_z[pointIndexOne] + o1->points_z[pointIndexTwo]) / 2;

		double temp = loopCenter.magnitude_x - o1->points_x[pointIndexOne];
		temp = temp * temp;
		double temp2 = loopCenter.magnitude_y - o1->points_y[pointIndexOne];
		temp2 = temp2 * temp2;
		double temp3 = loopCenter.magnitude_z - o1->points_z[pointIndexOne];
		temp3 = temp3 * temp3;
		temp = temp + temp2 + temp3;
		float a = (float)sqrt(temp);
		double axialOffset = o1->points_z[pointIndexOne];//z coord of the present loop[ ] TODO add axial offset back in?
		axialOffset = axialOffset - z1;//[ ] TODO check sign! field computed is relative to axialOffset
		float r = (float)sqrt(x1 * x1 + y1 * y1);
		if (axialOffset < grid_step/32)//if the plane of the field point is close to the plane of the loop: skip this point
		{//[ ] TODO tighten up this condition
			return;//result is 0 so far
		}
		loopCenterToFieldPoint.magnitude_x = x1 - loopCenter.magnitude_x;
		loopCenterToFieldPoint.magnitude_y = y1 - loopCenter.magnitude_y;
		double radius_angle = atan(loopCenterToFieldPoint.magnitude_x / loopCenterToFieldPoint.magnitude_y);
		double alpha = r / a;
		double beta = axialOffset / a;
		double gamma = axialOffset / r;
		double al2 = alpha * alpha;
		double be2 = beta * beta;
		double q = pow(1 + alpha, 2) + be2;
		double rq = sqrt(q);
		double k = sqrt((alpha * 4) / q);
		double fk = completeEllipticIntegralFirstKind(k);
		double ek = completeEllipticIntegralSecondKind(k);
		double Ht = loopCurrent / (2.0 * a * PI * rq);
		double Hx = Ht * (ek * (1 - al2 - be2) / (q - (alpha * 4)) + fk);
		double Hr = (r == 0.0) ? (0.0) : (Ht * gamma * (ek * (1 + al2 + be2) / (q - (alpha * 4)) - fk));
		fieldAtPointBxAxial = Hx * MU_0;//convert from H to B field
		fieldAtPointBrRadial = Hr * MU_0;
		fieldAtPointXyz.magnitude_x = fieldAtPointBrRadial * sin(radius_angle);//[ ] TODO check on sin, cos
		fieldAtPointXyz.magnitude_y = fieldAtPointBrRadial * cos(radius_angle);//make sure not backwards
		fieldAtPointXyz.magnitude_z = fieldAtPointBxAxial;
		resultantFieldXyz.magnitude_x = resultantFieldXyz.magnitude_x + fieldAtPointXyz.magnitude_x;
		resultantFieldXyz.magnitude_y = resultantFieldXyz.magnitude_y + fieldAtPointXyz.magnitude_y;
		resultantFieldXyz.magnitude_z = resultantFieldXyz.magnitude_z + fieldAtPointXyz.magnitude_z;

	}
	//2nd quark's current loops
	for (int i = 0; i < num_points_sphere; i += num_points_circle)
	{
		pointIndexOne = i;
		pointIndexTwo = pointIndexOne + num_points_circle / 2;
		loopCenter.magnitude_x = (o2->points_x[pointIndexOne] + o2->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o2->points_y[pointIndexOne] + o2->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o2->points_z[pointIndexOne] + o2->points_z[pointIndexTwo]) / 2;
		loopCharge = 0.0;
		for (int j = 0; j < num_points_circle; j++)
		{
			loopCharge += o2->charges[i + j];
		}
		loopCurrent = loopCharge * revsPerSecond; //Coulombs per second (Amps)

	//The following code is from main() of ELLIPTIC.C.  There are modifications, but the unmodified code is at the end of this file in a comment.
	//	/* input parameters */
		loopCenter.magnitude_x = (o2->points_x[pointIndexOne] + o2->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o2->points_y[pointIndexOne] + o2->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o2->points_z[pointIndexOne] + o2->points_z[pointIndexTwo]) / 2;

		double temp = loopCenter.magnitude_x - o2->points_x[pointIndexOne];
		temp = temp * temp;
		double temp2 = loopCenter.magnitude_y - o2->points_y[pointIndexOne];
		temp2 = temp2 * temp2;
		double temp3 = loopCenter.magnitude_z - o2->points_z[pointIndexOne];
		temp3 = temp3 * temp3;
		temp = temp + temp2 + temp3;
		loopCenterToFieldPoint.magnitude_x = x1 - loopCenter.magnitude_x;
		loopCenterToFieldPoint.magnitude_y = y1 - loopCenter.magnitude_y;
		double radius_angle = atan(loopCenterToFieldPoint.magnitude_x / loopCenterToFieldPoint.magnitude_y);
		float a = (float)sqrt(temp);
		double axialOffset = o2->points_z[pointIndexOne];//z coord of the present loop[ ] TODO add axial offset back in?
		axialOffset = axialOffset - z1;//[ ] TODO check sign! field computed is relative to axialOffset
		float r = (float)sqrt(x1 * x1 + y1 * y1);
		double alpha = r / a;
		double beta = axialOffset / a;
		double gamma = axialOffset / r;
		double al2 = alpha * alpha;
		double be2 = beta * beta;
		double q = pow(1 + alpha, 2) + be2;
		double rq = sqrt(q);
		double k = sqrt((alpha * 4) / q);
		double fk = completeEllipticIntegralFirstKind(k);
		double ek = completeEllipticIntegralSecondKind(k);
		double Ht = loopCurrent / (2.0 * a * PI * rq);
		double Hx = Ht * (ek * (1 - al2 - be2) / (q - (alpha * 4)) + fk);
		double Hr = (r == 0.0) ? (0.0) : (Ht * gamma * (ek * (1 + al2 + be2) / (q - (alpha * 4)) - fk));
		fieldAtPointBxAxial = Hx * MU_0;//convert from H to B field
		fieldAtPointBrRadial = Hr * MU_0;
		fieldAtPointXyz.magnitude_x = fieldAtPointBrRadial * sin(radius_angle);//[ ] TODO check on sin, cos
		fieldAtPointXyz.magnitude_y = fieldAtPointBrRadial * cos(radius_angle);//make sure not backwards
		fieldAtPointXyz.magnitude_z = fieldAtPointBxAxial;
		resultantFieldXyz.magnitude_x = resultantFieldXyz.magnitude_x + fieldAtPointXyz.magnitude_x;
		resultantFieldXyz.magnitude_y = resultantFieldXyz.magnitude_y + fieldAtPointXyz.magnitude_y;
		resultantFieldXyz.magnitude_z = resultantFieldXyz.magnitude_z + fieldAtPointXyz.magnitude_z;

	}

	//3rd quark's current loops
	for (int i = 0; i < num_points_sphere; i += num_points_circle)
	{
		pointIndexOne = i;
		pointIndexTwo = pointIndexOne + num_points_circle / 2;
		loopCenter.magnitude_x = (o3->points_x[pointIndexOne] + o3->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o3->points_y[pointIndexOne] + o3->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o3->points_z[pointIndexOne] + o3->points_z[pointIndexTwo]) / 2;
		loopCharge = 0.0;
		for (int j = 0; j < num_points_circle; j++)
		{
			loopCharge += o3->charges[i + j];
		}
		loopCurrent = loopCharge * revsPerSecond; //Coulombs per second (Amps)

	//The following code is from main() of ELLIPTIC.C.  There are modifications, but the unmodified code is at the end of this file in a comment.
	//	/* input parameters */
		loopCenter.magnitude_x = (o3->points_x[pointIndexOne] + o3->points_x[pointIndexTwo]) / 2;
		loopCenter.magnitude_y = (o3->points_y[pointIndexOne] + o3->points_y[pointIndexTwo]) / 2;
		loopCenter.magnitude_z = (o3->points_z[pointIndexOne] + o3->points_z[pointIndexTwo]) / 2;

		double temp = loopCenter.magnitude_x - o3->points_x[pointIndexOne];
		temp = temp * temp;
		double temp2 = loopCenter.magnitude_y - o3->points_y[pointIndexOne];
		temp2 = temp2 * temp2;
		double temp3 = loopCenter.magnitude_z - o3->points_z[pointIndexOne];
		temp3 = temp3 * temp3;
		temp = temp + temp2 + temp3;
		loopCenterToFieldPoint.magnitude_x = x1 - loopCenter.magnitude_x;
		loopCenterToFieldPoint.magnitude_y = y1 - loopCenter.magnitude_y;
		double radius_angle = atan(loopCenterToFieldPoint.magnitude_x / loopCenterToFieldPoint.magnitude_y);
		float a = (float)sqrt(temp);
		double axialOffset = o3->points_z[pointIndexOne];//z coord of the present loop[ ] TODO add axial offset back in?
		axialOffset = axialOffset - z1;//[ ] TODO check sign! field computed is relative to axialOffset
		float r = (float)sqrt(x1 * x1 + y1 * y1);
		double alpha = r / a;
		double beta = axialOffset / a;
		double gamma = axialOffset / r;
		double al2 = alpha * alpha;
		double be2 = beta * beta;
		double q = pow(1 + alpha, 2) + be2;
		double rq = sqrt(q);
		double k = sqrt((alpha * 4) / q);
		double fk = completeEllipticIntegralFirstKind(k);
		double ek = completeEllipticIntegralSecondKind(k);
		double Ht = (-1.0*loopCurrent) / (2.0 * a * PI * rq);//minus loop current for other direction of 2nd down quark
		double Hx = Ht * (ek * (1 - al2 - be2) / (q - (alpha * 4)) + fk);
		double Hr = (r == 0.0) ? (0.0) : (Ht * gamma * (ek * (1 + al2 + be2) / (q - (alpha * 4)) - fk));
		fieldAtPointBxAxial = Hx * MU_0;//convert from H to B field
		fieldAtPointBrRadial = Hr * MU_0;
		fieldAtPointXyz.magnitude_x = fieldAtPointBrRadial * sin(radius_angle);//[ ] TODO check on sin, cos
		fieldAtPointXyz.magnitude_y = fieldAtPointBrRadial * cos(radius_angle);//make sure not backwards
		fieldAtPointXyz.magnitude_z = fieldAtPointBxAxial;
		resultantFieldXyz.magnitude_x = resultantFieldXyz.magnitude_x + fieldAtPointXyz.magnitude_x;
		resultantFieldXyz.magnitude_y = resultantFieldXyz.magnitude_y + fieldAtPointXyz.magnitude_y;
		resultantFieldXyz.magnitude_z = resultantFieldXyz.magnitude_z + fieldAtPointXyz.magnitude_z;

	}
	//[X] TODO how to turn the down quark rotation in either direction? just flip current?
	//(did this with the term (-1.0*loopCurrent) in Ht above)

	result->magnitude_x = resultantFieldXyz.magnitude_x;
	result->magnitude_y = resultantFieldXyz.magnitude_y;
	result->magnitude_z = resultantFieldXyz.magnitude_z;
	result->magnitude = sqrt(result->magnitude_x * result->magnitude_x + result->magnitude_y * result->magnitude_y + result->magnitude_z * result->magnitude_z);
}

double completeEllipticIntegralFirstKind(double k)
{
	int ierr = 0;
	double firstKind = drf(0.0, 1.0 - pow(k, 2), 1.0, &ierr);
	if (ierr)
	{
		printf("invalid argument in drf\n");
	}
	return firstKind;
}
double completeEllipticIntegralSecondKind(double k)
{
	int ierr = 0;
	double secondKind = (drf(0.0, 1.0 - pow(k, 2), 1.0, &ierr)\
		- (pow(k, 2) / 3.0) * drd(0.0, 1.0 - pow(k, 2), 1.0, &ierr));
	if (ierr)
	{
		printf("invalid argument in drd\n");
	}
	return secondKind;


}


//---------------------------------------------------------------------------


/*
 *
 * ELLIPTIC.C:  Functions for computing complete and incomplete elliptic
 *              integrals of the first and second kind, and an example
 *              of usage for computing magnetic fields due to a circular
 *              current filament.
 *
 * DISCLAIMER
 *
 * Although every effort has been expended to ensure correctness, this
 * software is not guaranteed to be correct.
 *
 * You may do anything you want with this software, as long as you maintain
 * the existing comments intact, including the names of the original
 * authors.  If you modify this software, please indicate that fact in
 * your comments.
 *
 * E. Dennison (Nov 17, 1998)
 *
 */
#include "float.h"
#include "math.h"
#include "stdio.h"



#define MAX(x,y) (x>y?x:y)
#define MAX3(x,y,z) MAX(MAX(x,y),z)
#define MIN(x,y) (x>y?y:x)
#define MIN3(x,y,z) MIN(MIN(x,y),z)


 /*
  * NOTE: constants declared in float.h for the x86 CPU are:
  *
  * DBL_EPSILON	2.2204460492503131e-016,smallest such that 1.0+DBL_EPSILON!=1.0
  * DBL_MAX 	1.7976931348623158e+308  maximum possible value for type double
  * DBL_MIN 	2.2250738585072014e-308  minimum possible positive value for double
  *
  * If you are compiling this code for a non-Intel CPU, your values for these
  * constants may be different.
  */



  /*
   * F(k): complete elliptic integral of the first kind
   */

#define F(k,ierr) drf(0.0,1.0-pow(k,2),1.0,&ierr)

   /*
	* drf.c:   Compute the complete or incomplete elliptic integral of the
	*          first kind.
	*
	* Description:
	*
	*  For x,y and z non-negative and at most one of them zero, drf(x,y,z)
	*  = integral from zero to infinity of
	*
	*
	*            -1/2     -1/2     -1/2
	*  (1/2)(t+x)    (t+y)    (t+z)    dt.
	*
	*  If x, y or z is zero, the integral is complete.
	*
	*  *piErr returns non-zero if any arguments are invalid.
	*
	* Credits:
	*
	*  This function is adapted by E. Dennison from FORTRAN code written by:
	*
	*  Carlson, B.C.
	*  Notis, E.M
	*  Pexton, R.L.
	*
	*/
double drf(double x, double y, double z, int* piErr)
{
	int    iErr = 0;
	double mu,
		xn, yn, zn,
		xndev, yndev, zndev,
		xnroot, ynroot, znroot,
		lambda,
		epsilon,
		e2, e3,
		result,
		s;

	const double c1 = 1.0 / 24.0;
	const double c2 = 3.0 / 44.0;
	const double c3 = 1.0 / 14.0;
	const double errtol = pow(DBL_EPSILON * 4.0, 1.0 / 6.0);
	const double lolim = 5.0 * DBL_MIN;
	const double hilim = DBL_MAX / 5.0;

	if (piErr)
	{
		if (MIN3(x, y, z) < 0.0)
		{
			iErr = 1;
		}
		else if (MIN3(x + y, x + z, y + z) < lolim)
		{
			iErr = 2;
		}
		else if (MAX3(x, y, z) > hilim)
		{
			iErr = 3;
		}
	}
	if (iErr)
	{
		if (piErr)
		{
			*piErr = iErr;
		}
		result = 0.0;
	}
	else
	{
		xn = x;
		yn = y;
		zn = z;

		while (1)
		{
			mu = (xn + yn + zn) / 3.0;
			xndev = 2.0 - (mu + xn) / mu;
			yndev = 2.0 - (mu + yn) / mu;
			zndev = 2.0 - (mu + zn) / mu;
			epsilon = MAX3(fabs(xndev), fabs(yndev), fabs(zndev));
			if (epsilon < errtol) break;
			xnroot = sqrt(xn);
			ynroot = sqrt(yn);
			znroot = sqrt(zn);
			lambda = xnroot * (ynroot + znroot) + ynroot * znroot;
			xn = (xn + lambda) * 0.25;
			yn = (yn + lambda) * 0.25;
			zn = (zn + lambda) * 0.25;
		}
		e2 = xndev * yndev - pow(zndev, 2);
		e3 = xndev * yndev * zndev;
		s = 1.0 + (c1 * e2 - 0.1 - c2 * e3) * e2 + c3 * e3;

		if (piErr)
		{
			*piErr = 0;
		}
		result = s / sqrt(mu);
	}
	return result;
}
/* END drf() */

/*
 * E(k): complete elliptic integral of the second kind
 */

#define E(k,ierr)   (drf(0.0,1.0-pow(k,2),1.0,&ierr)\
                    -(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))

 /*
  * FastE(k): fast, complete elliptic integral of the second kind.
  *           Use this macro if the complete elliptic integral of
  *           the first kind was previously computed for the same
  *           value of k.
  */

#define FastE(F,k,ierr)	((F)-(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))

  /*
   * drd.c:   Compute the complete or incomplete elliptic integral of the
   *          second kind.
   *
   * Description:
   *
   *  For x and y non-negative, x+y and z positive, drf(x,y,z) = integral
   *  from zero to infinity of
   *
   *
   *            -1/2     -1/2     -3/2
   *  (3/2)(t+x)    (t+y)    (t+z)    dt.
   *
   *  If x or y is zero, the integral is complete.
   *
   *  *piErr returns non-zero if any arguments are invalid.
   *
   * Credits:
   *
   *  This function is adapted by E. Dennison from FORTRAN code written by:
   *
   *  Carlson, B.C.
   *  Notis, E.M
   *  Pexton, R.L.
   *
   */
double drd(double x, double y, double z, int* piErr)
{
	int     iErr = 0;
	double  mu,
		xn, yn, zn,
		xndev, yndev, zndev,
		xnroot, ynroot, znroot,
		lambda,
		epsilon,
		ea, eb, ec, ed, ef,
		sigma,
		power4,
		result,
		s1, s2;

	const double c1 = 3.0 / 14.0;
	const double c2 = 1.0 / 6.0;
	const double c3 = 9.0 / 22.0;
	const double c4 = 3.0 / 26.0;
	const double errtol = pow(DBL_EPSILON / 3.0, 1.0 / 6.0);
	double uplim;
	const double lolim = 2.0 / pow(DBL_MAX, 2.0 / 3.0);
	double tuplim = pow(DBL_MIN, 1.0 / 3.0);
	tuplim = pow(0.1 * errtol, 1.0 / 3.0) / tuplim;
	uplim = pow(tuplim, 2.0);

	if (piErr)
	{
		if (MIN(x, y) < 0.0)
		{
			iErr = 1;
		}
		else if (MAX3(x, y, z) > uplim)
		{
			iErr = 2;
		}
		else if (MIN(x + y, z) < lolim)
		{
			iErr = 3;
		}
	}
	if (iErr)
	{
		if (piErr)
		{
			*piErr = iErr;
		}
		result = 0.0;
	}
	else
	{
		xn = x;
		yn = y;
		zn = z;
		sigma = 0.0;
		power4 = 1.0;
		while (1)
		{
			mu = (xn + yn + 3.0 * zn) * 0.2;
			xndev = (mu - xn) / mu;
			yndev = (mu - yn) / mu;
			zndev = (mu - zn) / mu;
			epsilon = MAX3(fabs(xndev), fabs(yndev), fabs(zndev));
			if (epsilon < errtol) break;
			xnroot = sqrt(xn);
			ynroot = sqrt(yn);
			znroot = sqrt(zn);
			lambda = xnroot * (ynroot + znroot) + ynroot * znroot;
			sigma = sigma + power4 / (znroot * (zn + lambda));
			power4 = power4 * 0.25;
			xn = (xn + lambda) * 0.25;
			yn = (yn + lambda) * 0.25;
			zn = (zn + lambda) * 0.25;
		}
		ea = xndev * yndev;
		eb = zndev * zndev;
		ec = ea - eb;
		ed = ea - 6.0 * eb;
		ef = ed + ec + ec;
		s1 = ed * (-c1 + 0.25 * c3 * ed - 1.5 * c4 * zndev * ef);
		s2 = zndev * (c2 * ef + zndev * (-c3 * ec + zndev * c4 * ea));
		if (piErr)
		{
			*piErr = 0;
		}
		result = 3.0 * sigma + power4 * (1.0 + s1 + s2) / (mu * sqrt(mu));
	}
	return result;
}
/* END drd() */


/*
 * MAIN:    Sample program using drd and drf
 *
 * Description:
 *
 *          Demonstrate use of drd and drf in computing the magnetic
 *          field at any point in space due to a circular current
 *          filament (one turn coil).
 *
 * Author:
 *          E Dennison
 *
 * Disclaimer:
 *
 */
//int main()
//{
//#define PI (3.141592654)
//#define MU0 (PI*4.0E-7)
//
//	/* input parameters */
//	float a; /* loop radius */
//	float r; /* measurement point radius */
//	float x; /* measurement point axial position */
//	float i; /* loop current */
//	/* output parameters */
//	double Hx, Hr; /* axial, radial field components, A/m */
//	/* working vars */
//	double k, q, rq, fk, ek, al, be, ga, al2, be2, alt4, Ht;
//	int    ierr;
//
//	/* gather input parameters */
//	printf("Loop radius (meters): ");
//	scanf("%f", &a);
//	printf("Radius at measurement point (meters): ");
//	scanf("%f", &r);
//	printf("Axial position at measurement point (meters): ");
//	scanf("%f", &x);
//	printf("Loop current (amperes): ");
//	scanf("%f", &i);
//	/* begin computation here */
//	al = r / a;
//	alt4 = al * 4.0;
//	be = x / a;
//	ga = x / r;
//	al2 = al * al;
//	be2 = be * be;
//	q = pow(1 + al, 2) + be2;
//	rq = sqrt(q);
//	k = sqrt(alt4 / q);
//	fk = F(k, ierr);
//	ek = FastE(fk, k, ierr);
//	Ht = i / (2.0 * a * PI * rq);
//	Hx = Ht * (ek * (1 - al2 - be2) / (q - alt4) + fk);
//	Hr = (r == 0.0) ? (0.0) : (Ht * ga * (ek * (1 + al2 + be2) / (q - alt4) - fk));
//	/* display the output */
//	printf("Axial field,  Bx: %e (Tesla)\n", (float)Hx * MU0);
//	printf("Radial field, Br: %e (Tesla)\n", (float)Hr * MU0);
//
//	return 0;
//}
