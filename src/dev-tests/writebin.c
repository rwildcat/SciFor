#include<stdlib.h>
#include<stdio.h>

int main() {

	float x[10][20], nc;
	int i,j;
	FILE *fp;
	
	for(i=0; i<10; i++)
		for(j=0; j<20; j++)
			x[i][j] = j;
	
	nc = 21.;
	printf("n cols =%g\n", nc);
	
	fp = fopen("data.bin","wb");  // r for read, b for binary
	
	// write ndata
	fwrite( &nc, sizeof(nc), 1, fp); // write nc

	// write data
	fwrite( x, sizeof(x), 1, fp); // write x

	fclose(fp);
	
	return 0;

}