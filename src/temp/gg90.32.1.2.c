#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#define size 32
//#define PI 3.141592654
//#define SQH 0.707106
#define b 4
double temp1[size];
double pout[size];

inline sumdiff (double xy[size], int r);
 
main (int argc, char *argv[])
{
	static double SQH = M_SQRT2 / 2;
	
	int myrank, tag = 1234, count, numprocs;
	double mytime;
	MPI_Status status;
	int n=0, i, j, cj, m, h, k, sk, o1, vf;
	double angnum, ang, mult, Cos[size], Sin[size], temp[size], v[size],
    p1[size], p2[size], p1new[size], p3a1[size],p2b1[size],p2b2[size];
	double c, s, p1a[size], p2a[size], xnew1[size], xtemp1[size], p11[size],
    p12[size], pxnew[size], p2a1[size],p1a1[size],p1a2[size],x1[size],p2a2[size];
	double x[size], p21[size], p22[size], p1b[size], p2b[size], p2new[size], pd[size];
	double final[size], cfin[size],pnew1[size],pnew2[size],p1b1[size],p1b2[size];
	
	
	
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);

	// moved this calc of Cos and Sin tables to before the start time
	for (j = 1; j <= (n / 4); j++) {
		angnum = 4 * (j - 1) + 1;
		mult = M_PI / (n * 4);
		ang = angnum * mult;
		Cos[j] = cos (ang);
		Sin[j] = sin (ang);
	}
	
	// save the cost of doing this division through the program; store as a static for speed
	static int halfsize = size / 2;

	
	// move gg90 precalc here
	int m2 = 8;
	int j2 = 1, mj2=0;
	double Cm2[size], Sm2[size], angles2;
	for (j2 = 1; j2 <= (m2 / 4); j2++)
	{
		mj2 = 4 * (j2 - 1) + 1;
		angles2 = (M_PI * mj2) / (2 * m2);
		Cm2[j2] = cos (angles2);
		Sm2[j2] = sin (angles2);
	}

	
	mytime = MPI_Wtime ();
	
	if (myrank == 0) //procs 0 starts in here
	{
		//printf ("\n N = %d case on a Single processor\n", size);
		
		do {
			x[n] =n+1;
			// printf (" %lf\n", x[n]);
			n++;
		} while(n<size);
		///printf("n value is %d\n",n);
		
		for (j = 1; j <= (n / 4); j++)
		{
			// the following stuff was moved up before the timing
			// the assignment to c and s has been eliminated
//			angnum = 4 * (j - 1) + 1;
//			mult = PI / (n * 4);
//			ang = angnum * mult;
//			Cos[j] = cos (ang);
//			Sin[j] = sin (ang);
//			c = Cos[j];
//			s = Sin[j];
			
			cj = 2 * (j - 1);
			m = n - 1;
			
			temp[0] = x[cj];
			temp[1] = x[(m - cj)];
			temp[2] = x[(cj + (n / 2))];
			temp[3] = x[(m - (cj + (n / 2)))];
			
			//        mj = 4 * (j - 1) + 1;
			v[0] = ((Cos[j] * temp[0]) + (Sin[j] * temp[1]));
			v[1] = ((Sin[j] * temp[0]) - (Cos[j] * temp[1]));
			
			temp[5] = (temp[2] + temp[3]) * SQH;
			temp[6] = (temp[2] - temp[3]) * SQH;
			
			v[2] = ((Cos[j] * temp[5]) - (Sin[j] * temp[6]));
			v[3] = ((Sin[j] * temp[5]) + (Cos[j] * temp[6]));
			
			
			
			sumdiff (v, b);
			p1[(2 * j) - 2] = temp1[0];
			p1[(2 * j) - 1] = temp1[1];
			p2[(2 * j) - 2] = temp1[2];
			p2[(2 * j) - 1] = temp1[3];
		}
		
		//printf("pnew1 array is");
		
		sumdiff (p1, halfsize);	//n value by 2 half of the array
		for(i=0;i<8;i++)
			p11[i] = temp1[i];
		
		
		
		for(k=0,i=8;i<16;i++,k++)
			p12[k] = temp1[i];
		
		
		sumdiff (p11, k);
		
		for(i=0;i<8;i++)
			p1a[i] = temp1[i];//first 8 values are here
		
		for(i=0,k=4;k<8;k++,i++)
			pnew1[i] = p1a[k];
		
		gg90(pnew1,4, &Cm2, &Sm2);
		
		for(i=0,k=4;k<8;k++,i++) // remaining 4 values are here
			p1a[k] = pout[i];
		//p1a has 8 values perfectly in that
		//printf("p1a values are here\n");
		//for(i=0;i<8;i++)
		// printf("%lf\n",p1a[i]);
		
		sumdiff(p1a,b);
		p1a1[0] = temp1[0];
		p1a1[1] = temp1[1];
		p1a1[2] = (temp1[2] + temp1[3]) * SQH;
		p1a1[3] = (temp1[2] - temp1[3]) * SQH;
		//printf("\np1a2 values are here\n");
		//for(i=0;i<4;i++)
		// printf("%lf\n",p1a1[i]);
		
		for(i=0,k=4;k<8;k++,i++)
			p1a2[i] = p1a[k]; //second half
		
		sumdiff(p1a2,b);
		
		p1a2[0] = temp1[0];
		p1a2[1] = temp1[1];
		p1a2[2] = (temp1[2] + temp1[3]) * SQH;
		p1a2[3] = (temp1[2] - temp1[3]) * SQH;
		
		//for(i=0;i<4;i++)
		// printf("%lf\n",p1a2[i]);
		
		//assigning values in this way
		
		x1[0]  = p1a1[0];
		x1[31] = p1a1[1];
		x1[15] = p1a1[2];
		x1[16] = p1a1[3];
		
		x1[7]  = p1a2[0];
		x1[24]  = p1a2[1];
		x1[8]  = p1a2[2];
		x1[23]  = p1a2[3];
		
		gg90(p12,8, &Cm2, &Sm2);//need to make this as gg902
		
		
		sumdiff(pout,8);
		//printf("p1b\n\n");
		
		
		for(k=4,i=0;k<8;k++,i++)
		{
			p1b[i] = temp1[k];//4,7 values are here
			// printf("%lf\n",p1b[i]);
		}/*
		 gg90(p1b,4); //p1b(5:8)
		 
		 for(k=4,i=0;i<4;i++)
		 p1b[k] = pout[i]; //catch u ater
		 */
		
		for (k = 0; k < 4; k++)
			p1b1[k] = temp1[k];
		
		sumdiff(p1b1,b);
		
		p1b1[0] = temp1[0];
		p1b1[1] = temp1[1];
		p1b1[2] = (temp1[2] + temp1[3]) * SQH;
		p1b1[3] = (temp1[2] - temp1[3]) * SQH;
		
		//printf("\n\n p1b1 value is\n ");
		//for(i=0;i<4;i++)
		// printf("%lf\n",p1b1[i]);
		
		x1[3] = p1b1[0];
		x1[28] = p1b1[1];
		x1[12] = p1b1[2];
		x1[19] = p1b1[3];
		
		gg90(p1b,4,&Cm2, &Sm2); //p1b(5:8)
		
		for(i=0;i<4;i++)
			p1b[i] = pout[i]; //catch u ater
		
		sumdiff(p1b,b);
		
		p1b2[0] = temp1[0];
		p1b2[1] = temp1[1];
		p1b2[2] = (temp1[2] + temp1[3]) * SQH;
		p1b2[3] = (temp1[2] - temp1[3]) * SQH;
		
		//printf("\np1b2 value is\n ");
		//for(i=0;i<4;i++)
		//  printf("%lf\n",p1b2[i]);
		
		x1[4] = p1b2[0];
		x1[27] = p1b2[1];
		x1[11] = p1b2[2];
		x1[20] = p1b2[3];
		
		//printf("\n****************\n");
		//printf("X values are\n");
		//for(i=0;i<32;i++)
		// printf("%lf\n",x1[i]);
		
		gg90(p2,halfsize, &Cm2, &Sm2);
		
		sumdiff(pout, halfsize);
		
		for(i=0;i<8;i++)
			p21[i] = temp1[i];
		
		for(i=0,k=8;k<16;k++,i++)
			p22[i] = temp1[k];
		
		sumdiff(p21,8);
		
		//printf("processor2a\n");
		//  for (i = 0; i < 8; i++)
		//  printf("%lf\n",temp1[i]);
		
		for(k=4,i=0;k<8;k++,i++)
			p2a[i] = temp1[k];//4,7 values are here
		
		for (k = 0; k < 4; k++)
			p2a1[k] = temp1[k];
		
		sumdiff(p2a1,b);
		
		p2a1[0] = temp1[0];
		p2a1[1] = temp1[1];
		p2a1[2] = (temp1[2] + temp1[3]) * SQH;
		p2a1[3] = (temp1[2] - temp1[3]) * SQH;
		
		
		gg90(p2a,4, &Cm2, &Sm2); //p1b(5:8)
		
		sumdiff(pout,4);
		
		p2a2[0] = temp1[0];
		p2a2[1] = temp1[1];
		p2a2[2] = (temp1[2] + temp1[3]) * SQH;
		p2a2[3] = (temp1[2] - temp1[3]) * SQH;
		
		x1[1]  = p2a1[0];
		x1[30] = p2a1[1];
		x1[14] = p2a1[2];
		x1[17] = p2a1[3];
		
		x1[6]  = p2a2[0];
		x1[25] = p2a2[1];
		x1[9]  = p2a2[2];
		x1[22] = p2a2[3];
		
		//printf("\n****************\n");
		
		gg90(p22,8, &Cm2, &Sm2);
		sumdiff(pout,8);
		
		for(k=4,i=0;k<8;k++,i++)
			p2b[i] = temp1[k];//4,7 values are here
		
		
		for (k = 0; k < 4; k++)
			p2b1[k] = temp1[k];
		
		sumdiff(p2b1,b);
		
		p2b1[0] = temp1[0];
		p2b1[1] = temp1[1];
		p2b1[2] = (temp1[2] + temp1[3]) * SQH;
		p2b1[3] = (temp1[2] - temp1[3]) * SQH;
		
		gg90(p2b,b, &Cm2, &Sm2);
		
		sumdiff(pout,4);
		
		p2b2[0] = temp1[0];
		p2b2[1] = temp1[1];
		p2b2[2] = (temp1[2] + temp1[3]) * SQH;
		p2b2[3] = (temp1[2] - temp1[3]) * SQH;
		
		x1[2]  = p2b1[0];
		x1[29] = p2b1[1];
		x1[13] = p2b1[2];
		x1[18] = p2b1[3];
		
		x1[5]  = p2b2[0];
		x1[26] = p2b2[1];
		x1[10] = p2b2[2];
		x1[21] = p2b2[3];
		
		//     printf("\n****************\n");
		//    printf("X values are\n");
		//   for(i=0;i<32;i++)
		// printf("%lf\n",x1[i]);
		
		mytime = MPI_Wtime () - mytime;	//get the time just after work is done and take the difference 
		mytime = mytime * 1000000;
		printf ("Timing from rank %d is %lfus.\n", myrank, mytime);
	}//end of rank 2
	MPI_Finalize ();
	return 0;
 	
}	//end of main program			


// make the sumdiff inline to ensure it appears inline, if possible
// merge four loops into two
// tried making loops countdown to a test !=0 but it didn't seem to help
inline sumdiff (double xy[size], int r)
{
	int f1, f2, i, k; // j, h, l
	double x1[size], x2[size];
	f1 = (r / 2);
//	f2 = (r / 2) + 1;

	
	for (i = 0; i < f1; i++) {
		x1[i] = xy[i];
		x2[i] = xy[i + f1];
	}

//	for (h = 0, j = f1; j < r; j++, h++)
//		x2[h] = xy[j];
	
	for (k = 0; k < f1; k++) {
		temp1[k] = x1[k] + x2[k];
		temp1[k+f1] = x1[k] - x2[k];
	}
	
//	for (l = f1, k = 0; l < r; l++, k++)
//		temp1[l] = x1[k] - x2[k];

}

	
gg90 (double g2[size], int m2, double *Cm2, double *Sm2)
{
	int j2 = 1, mj2, i, p;	//m2 conatins teh length of array and g2 has the array elements of p2
	double c, s, angles2; // eliminate Cm2[size], Sm2[size]
	
	for (j2 = 1; j2 <= (m2 / 4); j2++)
    {
		mj2 = 4 * (j2 - 1) + 1;
		/* the following is now precalc
		angles2 = (M_PI * mj2) / (2 * m2);
		Cm2[j2] = cos (angles2);
		Sm2[j2] = sin (angles2);
		*/
		
		//c = Cm2[j2];
		//s = Sm2[j2];
		
		p = mj2; //4 * (j2 - 1) + 1;
		
		pout[mj2 - 1] = ((Cm2[j2] * g2[p - 1]) + (Sm2[j2] * g2[p]));	//overwriting the values check that
		pout[mj2] = (Sm2[j2] * g2[p - 1]) - (Cm2[j2] * g2[p]);
		pout[mj2 + 1] = (-Sm2[j2] * g2[p + 1]) + (Cm2[j2] * g2[p + 2]);
		pout[mj2 + 2] = (Cm2[j2] * g2[p + 1]) + (Sm2[j2] * g2[p + 2]);
		
    }
}

// orignal output: as is
// mod1 output: basic output with varying value for size, and improvements except c and s subst
// mod2 output: above, plus c and s subts; before inline
// mod3 output: above plus sumdiff inline
// mod4 ouput: replace PI and SQH
// mod5 output: precalc of gg90
// mod6 output: elim c and s assignment in gg90; p and mj2 have same assignment so only do once
// mod7 output; realized size was incorrectly set to 16, so modified to 32 and got drastic improvement
