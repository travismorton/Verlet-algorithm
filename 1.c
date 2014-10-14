#include <stdio.h>
#include <stdlib.h>
#include <time.h>

extern float pos[2][16];
extern float vel[2][16];

int main (){
	float pos[2][16] = {{0, 0, 0, 0,
			1, 1, 1, 1,
			2, 2, 2, 2,
			3, 3, 3, 3},
		    {0, 1, 2, 3,
			0, 1, 2, 3,
			0, 1, 2, 3,
			0, 1, 2, 3}};
	float vel[2][16];
	int i, j;
	
	srand48(1);
	for (i = 0; i < 2; i++){
		for (j = 0; j < 16; j++){
			vel[i][j] = drand48()*3.0-1.5;
			//printf("vel[%d][%d] = %f\n", i, j, vel[i][j]);
		}
	}
	return 0;
}
