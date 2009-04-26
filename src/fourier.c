#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#define PI 3.141592654
#define SQH 0.707106
int ordb[65536];
double sum2[65536];
double poutf[65536];

main ( int argc, char *argv[] )
{
    int myrank, numprocs,tag=1234,c1;
    int n, i, n1 = 0, m, space, nstart,nstart2,space2, ncounter, m2, j, k, h1;
    n = atoi ( argv[1] );
    double xout[n], J, mytime, x[n], x2[n], x4[n], temp[n], temp1[n], temp2[n], xbothalf[n], xbothalf2[n], temp5[n], temp6[n], x1out[n], x2out[n], fxout1[n],fxout2[n],fxout3[n],fxout4[n];

    void dct ( double [], int );
    int getlog ( int );

    MPI_Status status;
    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size ( MPI_COMM_WORLD, &numprocs );

    if ( myrank == 0 ) {
        fftporder4 ( n );

        while ( n1<n ) {
            x[n1] = n1+1;
            n1++;
        }
        mytime = MPI_Wtime ();

        J = getlog ( n );

        m = n;
        space = 2;
        nstart = 2;

        for ( ncounter = 1; ncounter <= ( J - 2 ); ncounter++ ) {
            sumdiff2 ( x, 0, m, n );

            for ( i=0;i<n;i++ )
                x[i]=sum2[i];

            m2 = m / 2;

            for ( i = m2, j = 1; i < m; i++, j++ )
                temp[j] = x[i];

            for ( k = 0, i = m2, j = 1; i >= ( m2 / 2 ) + 1; i--, j++, k++ )
                temp1[k] = temp[i] - temp[j];


            k = 2* ( n/4 );
            MPI_Send ( temp1, k, MPI_DOUBLE, 2, tag, MPI_COMM_WORLD );


            for ( k = 0, i = m2 / 2, j = ( m2 / 2 ) + 1; i >= 1, j <= m2;
                    i--, j++, k++ )
                temp2[k] = temp[i] + temp[j];

            k = 2* ( n/4 );
            MPI_Send ( temp2, k , MPI_DOUBLE, 3, tag, MPI_COMM_WORLD );

            //DCT IMPLEMENTATION

            MPI_Recv ( temp1, k , MPI_DOUBLE, 2, tag,MPI_COMM_WORLD, &status );

            MPI_Recv ( temp2, k , MPI_DOUBLE, 3, tag,MPI_COMM_WORLD, &status );

            for ( i = 0; i < ( m2 / 2 ); i++ ) {
                xbothalf[i] = temp2[i];
                xbothalf2[i] = -1 * temp1[i];
            }

            for ( i = ( m2 / 2 ), j = ( m2/2 )-1; i<m2; i++, j-- ) {
                xbothalf[i] = temp2[j];
                xbothalf2[i] = temp1[j];
            }

            nstart2 = nstart;
            space2 = space;

            for ( i = 0; i < m2,nstart<=n; i++ ) {
                ordb[i] = ( ordb[i] / 2 );
                fxout1[nstart] = xbothalf[ordb[i]];
                fxout2[nstart] = xbothalf2[ordb[i]];
                nstart=nstart+space;
            }

            nstart = nstart2 + ( space2 / 2 );
            m = m / 2;
            space = 2 * space2;
        }

        //Last Step 4 Values

        for ( i = 0; i < 4; i++ )
            temp1[i] = x[i];

        sumdiff2 ( temp1, 0, 4, n );

        for ( i = 0; i < 4; i++ )
            x[i] = sum2[i];

        sumdiff2 ( x, 0, 2, n );

        fxout1[1] = sum2[0];
        fxout2[1] =0;

        fxout1[1 + ( n / 2 ) ] = sum2[1];
        fxout2[1+ ( n/2 ) ]=0;

        for ( j = 0, i = 3; i <= 4; i++, j++ )
            temp[j] = x[i-1];

        fxout1[ ( n / 4 ) + 1] = temp[0];
        fxout2[ ( n / 4 ) + 1] = -1 * temp[1];

        fxout1[ ( n / 4 ) + ( n / 2 ) + 1] = temp[0];
        fxout2[ ( n / 4 ) + ( n / 2 ) + 1] = temp[1];

        MPI_Recv ( fxout3, n+1 , MPI_DOUBLE, 1, tag,MPI_COMM_WORLD, &status );
        MPI_Recv ( fxout4, n+1 , MPI_DOUBLE, 1, tag,MPI_COMM_WORLD, &status );

        for ( i = 1; i <= n; i++ ) {
            fxout1[i] = fxout1[i] - fxout4[i];
            fxout2[i] = fxout2[i] + fxout3[i];
        }

        //     for(i=1;i<=n;i++)
        //  printf("%lf %lfi\n",fxout1[i],fxout2[i]);

        mytime = MPI_Wtime () - mytime;
        mytime = mytime * 1000000;
        printf ( "\nTiming from rank %d is %lfus.\n", myrank,mytime );
    }   //if rank 0 ends here

    if ( myrank == 2 ) {
        k = ( n/4 ) *2;
        J = getlog ( n );
        J = J-2;
        for ( ncounter=J;ncounter>=1;ncounter-- ) {
            MPI_Recv ( temp1, k , MPI_DOUBLE, 0, tag,MPI_COMM_WORLD, &status );
            c1 =  pow ( 2, ncounter );
            dct ( temp1,c1 );
            MPI_Send ( temp1, k , MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );

        }
    }

    if ( myrank == 3 ) {
        k = ( n/4 ) *2;
        J = getlog ( n );
        J = J-2;

        for ( ncounter=J;ncounter>=1;ncounter-- ) {
            MPI_Recv ( temp2, k , MPI_DOUBLE, 0, tag,MPI_COMM_WORLD, &status );
            c1 =  pow ( 2, ncounter );
            dct ( temp2, c1 );
            MPI_Send ( temp2, k , MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );
        }
    }

    if ( myrank == 1 ) {
        fftporder4 ( n );

        while ( n1<n ) {
            x2[n1] = n1+1;
            n1++;
        }

        J = getlog ( n );

        m = n;
        space = 2;
        nstart = 2;

        for ( ncounter = 1; ncounter <= ( J - 2 ); ncounter++ ) {
            sumdiff2 ( x2, 0, m, n );

            for ( i=0;i<n;i++ )
                x2[i]=sum2[i];

            m2 = m / 2;

            for ( i = m2, j = 1; i < m; i++, j++ )
                temp[j] = x2[i];

            for ( k = 0, i = m2, j = 1; i >= ( m2 / 2 ) + 1; i--, j++, k++ )
                temp1[k] = temp[i] - temp[j];

            k = 2* ( n/4 );
            MPI_Send ( temp1, k, MPI_DOUBLE, 4, tag, MPI_COMM_WORLD );


            for ( k = 0, i = m2 / 2, j = ( m2 / 2 ) + 1; i >= 1, j <= m2;
                    i--, j++, k++ )
                temp2[k] = temp[i] + temp[j];

            k = 2* ( n/4 );
            MPI_Send ( temp2, k, MPI_DOUBLE, 5, tag, MPI_COMM_WORLD );

            MPI_Recv ( temp1, k , MPI_DOUBLE, 4, tag,MPI_COMM_WORLD, &status );

            MPI_Recv ( temp2, k , MPI_DOUBLE, 5, tag,MPI_COMM_WORLD, &status );

            for ( i = 0; i < ( m2 / 2 ); i++ ) {
                xbothalf[i] = temp2[i];
                xbothalf2[i] = -1 * temp1[i];
            }

            for ( i = ( m2 / 2 ), j = ( m2/2 )-1; i<m2; i++, j-- ) {
                xbothalf[i] = temp2[j];
                xbothalf2[i] = temp1[j];
            }

            nstart2 = nstart;
            space2 = space;

            for ( i = 0; i < m2,nstart<=n; i++ ) {
                ordb[i] = ( ordb[i] / 2 );
                fxout3[nstart] = xbothalf[ordb[i]];
                fxout4[nstart] = xbothalf2[ordb[i]];
                nstart=nstart+space;
            }

            nstart = nstart2 + ( space2 / 2 );
            m = m / 2;
            space = 2 * space2;
        }   //for counter endsa

        //Last Step 4 Values

        for ( i = 0; i < 4; i++ )
            temp1[i] = x2[i];

        sumdiff2 ( temp1, 0, 4, n );

        for ( i = 0; i < 4; i++ )
            x2[i] = sum2[i];

        sumdiff2 ( x2, 0, 2, n );

        fxout3[1] = sum2[0];
        fxout4[1] =0;

        fxout3[1 + ( n / 2 ) ] = sum2[1];
        fxout4[1+ ( n/2 ) ]=0;

        for ( j = 0, i = 3; i <= 4; i++, j++ )
            temp[j] = x2[i-1];

        fxout3[ ( n / 4 ) + 1] = temp[0];
        fxout4[ ( n / 4 ) + 1] = -1 * temp[1];

        fxout3[ ( n / 4 ) + ( n / 2 ) + 1] = temp[0];
        fxout4[ ( n / 4 ) + ( n / 2 ) + 1] = temp[1];

        MPI_Send ( fxout3, n+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );
        MPI_Send ( fxout4, n+1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );

    }    //if rank 1 ends here
    if ( myrank == 4 ) {
        k = ( n/4 ) *2;
        J = getlog ( n );
        J = J-2;
        for ( ncounter=J;ncounter>=1;ncounter-- ) {
            MPI_Recv ( temp1, k , MPI_DOUBLE, 1, tag,MPI_COMM_WORLD, &status );
            c1 =  pow ( 2, ncounter );
            dct ( temp1,c1 );
            MPI_Send ( temp1, k , MPI_DOUBLE, 1, tag, MPI_COMM_WORLD );

        }
    }

    if ( myrank == 5 ) {
        k = ( n/4 ) *2;
        J = getlog ( n );
        J = J-2;

        for ( ncounter=J;ncounter>=1;ncounter-- ) {
            MPI_Recv ( temp2, k , MPI_DOUBLE, 1, tag,MPI_COMM_WORLD, &status );
            c1 =  pow ( 2, ncounter );
            dct ( temp2, c1 );
            MPI_Send ( temp2, k , MPI_DOUBLE, 1, tag, MPI_COMM_WORLD );
        }
    }
    MPI_Finalize ();
    return 0;
}

//**************DCT IMPLEMENTATION*****************//
void dct ( double xout[], int z1 )
{
    long int n = 0;
    int L, m2, myrank, pcounter, l1, tag =1234, h1, h3, h4, h5, h6, K1, h16, h7, h8, numprocs, i, j, cj, m, k,
                                          count;
    long int j4r[z1];
    double c, s, J1, angnum, mrs, ang, mult, C[z1], S[z1], C2b[z1], S2b[z1],
    angle, angle1, temp[z1], xin[z1], temprs[z1], v[z1], S2[z1], p1[z1],
    C2c[z1], S2c[z1], temp2rs[z1], p2[z1], mytime, temp1[z1], block[z1],
    pair1[z1], pair2[z1], temp2[z1], temp3[z1], temp4[z1], p[z1],xpin[z1];

    while ( n < z1 ) {
        n++;
    }

    if ( n == 2 ) {
        angle = PI / 8;
        c = cos ( angle );
        s = sin ( angle );

        xpin[0] = xout[0];
        xpin[1] = xout[1];

        xout[0] = ( c * xpin[0] ) + ( s * xpin[1] );
        xout[1] = ( -s * xpin[0] ) + ( c * xpin[1] );
    } else if ( n == 4 ) {
        c = cos ( PI / 16 );
        s = sin ( PI / 16 );

        temp[1] = xout[0];
        temp[2] = xout[3];
        temp[3] = xout[2];
        temp[4] = xout[1];

        xin[0] = c * temp[1] + s * temp[2];
        xin[1] = ( -s * temp[1] ) + c * temp[2];

        c = cos ( ( 5 * PI ) / 16 );
        s = sin ( ( 5 * PI ) / 16 );

        xin[2] = c * temp[3] + s * temp[4];
        xin[3] = ( -s * temp[3] ) + c * temp[4];

        sumdiff2 ( xin, 0, 4, z1 );

        temp[0] = sum2[2] * SQH;
        temp[1] = sum2[3] * SQH;

        xout[0] = sum2[0];
        xout[1] = sum2[1];
        xout[2] = temp[0] + temp[1];
        xout[3] = -temp[0] + temp[1];
    } else {
        for ( j = 1; j <= ( n / 2 ); j++ ) {
            angnum = 4 * ( j - 1 ) + 1;
            mult = PI / ( n * 4 );
            ang = angnum * mult;
            C[j] = cos ( ang );
            S[j] = sin ( ang );
        }

        for ( j = 1; j <= ( n / 4 ); j++ ) {
            cj = 2 * ( j - 1 );
            m = n - 1;
            //0,7,4,3 in the first iteration and 2 5 6 1 in the second iteration
            pair1[0] = xout[cj];
            pair1[1] = xout[ ( m - cj ) ];
            pair2[0] = xout[ ( cj + ( n / 2 ) ) ]; //xtemp(3:4)
            pair2[1] = xout[ ( m - ( cj + ( n / 2 ) ) ) ];

            v[0] = ( C[j] * pair1[0] ) + ( S[j] * pair1[1] );
            v[1] = ( -S[j] * pair1[0] ) + ( C[j] * pair1[1] );

            v[2] = ( C[j + ( n / 4 ) ] * pair2[0] ) + ( S[j + ( n / 4 ) ] * pair2[1] );
            v[3] = ( -S[j + ( n / 4 ) ] * pair2[0] ) + ( C[j + ( n / 4 ) ] * pair2[1] );

            sumdiff2 ( v, 0, 4, z1 );

            p1[ ( 2 * j ) - 2] = sum2[0];
            p1[ ( 2 * j ) - 1] = sum2[1];
            p2[ ( 2 * j ) - 2] = sum2[2];
            p2[ ( 2 * j ) - 1] = sum2[3];
        }

        pcounter = 1;
        while ( pcounter <= 2 ) {
            if ( pcounter == 1 ) {
                for ( i = 0, j = 1; i < n / 2; i++, j++ ) {
                    p[i] = p1[i];
                    count++;
                }
            }   //if ends here
            else {
                count = 0;
                for ( j = 1, i = 0; i < n / 2; i++, j++ ) {
                    p[j] = p2[i];
                    count++;
                }
                for ( i = 1; ( 4 * ( i - 1 ) + 1 ) <= n / 2; i++ ) {
                    angle = 4 * ( i - 1 ) + 1;
                    angle1 = angle / n;
                    C2b[i] = cos ( PI * angle1 );
                    S2b[i] = sin ( PI * angle1 );
                }
                gg90cs ( p, C2b, S2b, count, z1 );

                for ( i = 0, j = 1; j <= n / 2; j++, i++ ) {
                    p[i] = poutf[j];
                }
            }   //else ends here

            J1 = getlog ( n );
            K1 = J1 - 2;
            m = n / 2;

            for ( l1 = 1; l1 <= K1; l1++ ) {
                m2 = pow ( 2, ( l1 - 1 ) );
                L = m / m2;

                for ( j = 1; ( 4 * ( j - 1 ) ) < L / 2; j++ ) {
                    mrs = 4 * ( j - 1 ) + 1;
                    mrs = mrs / L;
                    C2c[j] = cos ( PI * mrs );
                    S2c[j] = sin ( PI * mrs );
                }

                for ( k = 1; k <= m2; k++ ) {
                    h1 = L * ( k - 1 );
                    h3 = ( L * k ) - 1;
                    h6 = L * ( k - 1 );
                    h16 = 0;
                    while ( h1 <= h3 ) {
                        block[h1] = p[h1];
                        h16++;
                        h1++;
                    }
                    sumdiff2 ( block, h6, h3 + 1, z1 );

                    for ( i = 0; i < h16; i++ ) {
                        temp2[i] = sum2[i];
                    }

                    h4 = ( L / 2 );


                    for ( i = 1; h4 < L; h4++, i++ ) {
                        temp3[i] = temp2[h4];
                    }
                    gg90cs ( temp3, C2c, S2c, L / 2, z1 );

                    for ( i = 1, h5 = L / 2; h5 < L; i++, h5++ ) {
                        temp2[h5] = poutf[i];
                    }

                    h7 = L * ( k - 1 );
                    h8 = ( L * k ) - 1;

                    for ( i = 0; h7 <= h8; i++, h7++ ) {
                        p[h7] = temp2[i];
                    }
                }
                if ( pcounter == 1 ) {
                    for ( i = 0; i < n / 2; i++ )
                        xout[i] = p[i];
                } else {
                    for ( j = 0, i = n / 2; i < n; i++, j++ )
                        xout[i] = p[j];
                }
            }
            pcounter++;
        }

    }
}    //dct ends here

gg90cs ( double pin1[], double C2b[], double S2b[], int count, int z1 )
{
    int m, j, mj, mj2, i;
    double temp[z1], temp1[z1], pout3[z1], pout4[z1];
    m = count;

    if ( m == 2 ) {
        poutf[1] = ( pin1[1] + pin1[2] ) * S2b[1];
        poutf[2] = ( -pin1[1] + pin1[2] ) * S2b[1];
    } else {
        for ( j = 1; j <= m / 4; j++ ) {
            mj = ( 2 * j ) - 1;

            temp[1] = pin1[mj];
            temp[2] = pin1[mj + 1];

            pout3[mj] = ( C2b[j] * temp[1] ) + ( S2b[j] * temp[2] );
            pout3[mj + 1] = ( -S2b[j] * temp[1] ) + ( C2b[j] * temp[2] );

            temp1[1] = pin1[ ( ( m / 2 ) + mj ) ];
            temp1[2] = pin1[ ( ( m / 2 ) + mj + 1 ) ];

            pout4[mj] = ( -S2b[j] * temp1[1] ) + ( C2b[j] * temp1[2] );
            pout4[mj + 1] = ( -C2b[j] * temp1[1] ) + ( -S2b[j] * temp1[2] );
        }

        for ( i = 1; i <= m / 2; i++ ) {
            poutf[i] = pout3[i];
        }
        for ( i = 1, j = ( m / 2 + 1 ); i <= m; i++, j++ ) {
            poutf[j] = pout4[i];
        }
    }    //else ends
}

////////////////////////////////////
sumdiff2 ( double xy3[], int r1, int r2, int z1 )
{
    int f1, r, i, j, h, k, l, size1;
    size1 = z1;
    double x1[size1], x2[size1], xy2[size1];

    r = r2 - r1;
    f1 = ( r / 2 );

    for ( i = r1, k = 0; i < r2; i++, k++ )
        xy2[k] = xy3[i];

    for ( i = 0, k = 0; i < f1; i++, k++ )
        x1[k] = xy2[i];

    for ( h = 0, j = f1; j < r; j++, h++ )
        x2[h] = xy2[j];

    for ( k = 0; k < f1; k++ )
        sum2[k] = x1[k] + x2[k];

    for ( l = f1, k = 0; l < r; l++, k++ )
        sum2[l] = x1[k] - x2[k];
}
fftporder4 ( int n1 )
{
    double J;
    double veca[n1];
    veca[0] = 0;
    veca[1] = 2;
    veca[2] = 1;
    veca[3] = 3;
    double vecb[n1], z1[n1];
    int klev, i, nv, s1, k;
    nv = 2;
    int n=n1;
    J = getlog ( n );
    for ( klev = 1; klev <= J - 2; klev++ ) {
        nv *= 2;

        for ( i = 0; i <= ( n1 / 2 ) - 1; i++ ) {
            vecb[i] = 2 * veca[i];
        }

        s1 = ( 2 * nv ) - 1;

        for ( i = 0; i < nv; i++ )
            z1[i] = s1 - vecb[i];

        for ( i = nv, k = nv - 1; k >= 0, i < ( 2 * nv ); k--, i++ )
            vecb[i] = z1[k];

        for ( i = 0; i < ( nv * 2 ); i++ ) {
            veca[i] = vecb[i];
        }
    }

    for ( i = 0; i < n1; i++ ) {
        ordb[i] = vecb[i];
    }
}
int getlog ( int n )
{
    if ( n==8 ) return ( 3 );
    else if ( n==16 ) return ( 4 );
    else if ( n==32 ) return ( 5 );
    else if ( n==64 ) return ( 6 );
    else if ( n==128 ) return ( 7 );
    else if ( n==256 ) return ( 8 );
    else if ( n==512 ) return ( 9 );
    else if ( n==1024 ) return ( 10 );
    else if ( n==2048 ) return ( 11 );
    else if ( n==4096 ) return ( 12 );
    else if ( n==8192 ) return ( 13 );
    else if ( n==16384 ) return ( 14 );
    else if ( n==32768 ) return ( 15 );
    else if ( n==65536 ) return ( 16 );
    else  return ( 0 );
}