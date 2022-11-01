#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
//----------------------------------------------------------------------
#define nx  200
#define ny 	51
#define q		  9
//----------------------------------------------------------------------
const int time			=  100000;
const int noOfSnaps	= 		 2;
const int dispFreq	= 	  50;
const double tau		=		0.95;
const double rho0		= 	 1.0;
const double uInlet	=	 0.01;
//----------------------------------------------------------------------
int main(void)
{
	int i,j,a,a1,ts,ia,ja;
	static int ex[q],ey[q],kb[q];
	static int isn[nx+2][ny+2];
  double Q_xx[q],Q_yy[q],Q_xy[q],cc_xx[q],cc_xy[q],cc_yy[q];
	double tmp1, tmp2, tmp3, cs2Inv, rhoAvg, feq, f_neq_[q],f_neq, fx_t, fy_t, Cd, Cl;
	double dudx, dudy, dvdx, dvdy, cDOTu, uDOTu, term1, term2, term3;
	double drhou2dx, drhouvdx, drhouvdy, drhovvdy, dQddRhouudx, dQddRhouudy;
	double QddRhouu, QddRhouuxp1, QddRhouuxp2, QddRhouuym1, QddRhouuym2, QddRhouuyp1, QddRhouuyp2;
	double fx[2], fy[2];
	static double f[q][nx+2][ny+2], ft[q][nx+2][ny+2];
	static double wt[q], ux[nx+2][ny+2], uy[nx+2][ny+2], rho[nx+2][ny+2],ux_a[ny+2];

	char prefix[]="snap_",type[]=".dat",filename[15],solstr[5];
	int solnumber=0;
	cs2Inv=3.0;
	FILE *soln, *fid, *fid2, *fid3;

 for (j=1; j<=ny; j++)
 {
   ux_a[j]=6*uInlet*(ny-(j-0.5))*(j-0.5)/(ny*ny);
 }
	//----------------------------------------------------------------------
	ex[0] =	0; ey[0] = 0;
	ex[1] = 1; ey[1] = 0;
	ex[2] = 0; ey[2] = 1;
	ex[3] =-1; ey[3] = 0;
	ex[4] = 0; ey[4] =-1;
	ex[5] = 1; ey[5] = 1;
	ex[6] =-1; ey[6] = 1;
	ex[7] =-1; ey[7] =-1;
	ex[8] = 1; ey[8] =-1;
	//----------------------------------------------------------------------
	for (a=0; a<9; a++)
	{
		if (a==0) 				{wt[a] = 4.0/9.0 ;}
		if (a>=1 && a<=4)	{wt[a] = 1.0/9.0 ;}
		if (a>=5 && a<=8)	{wt[a] = 1.0/36.0;}
	}
	//----------------------------------------------------------------------
	for (a=0; a<q; a++)
	{
		for (a1=a; a1<q; a1++)
		{
			if ( ex[a]+ex[a1]==0 && ey[a]+ey[a1]==0)
			{
				kb[a]	= a1;
				kb[a1]=  a;
			}
		}
	}
	//----------------------------------------------------------------------
	for (i=0; i<=nx+1; i++)
	{
		for (j=0; j<=ny+1; j++)
		{
			for (a=0; a<9; a++)
			{
				cDOTu		= uInlet*ex[a];
				uDOTu		= uInlet*uInlet;
				f[a][i][j] = wt[a]*rho0;//*(1.0 + 3.0*cDOTu + 4.5*cDOTu*cDOTu - 1.5*uDOTu);
			}
		}
	}
	//----------------------------------------------------------------------
	fid = fopen("tRhoUx.dat","w");
	fprintf(fid,"Variables=time,rho,ux\n");
 	fid2 = fopen("xRhoUxUy.dat","w");
	fprintf(fid2,"Variables=x,rho,ux,uy\n");
 	fid3 = fopen("yUxaUx.dat","w");
	fprintf(fid3,"Variables=y,uxa,ux1,ux2,ux3\n");
	//----------------------------------------------------------------------
	for(a=0;a<9;a++)
	{
		Q_xx[a] = ex[a]*ex[a] - 1.0/cs2Inv;
		Q_yy[a] = ey[a]*ey[a] - 1.0/cs2Inv;
		Q_xy[a] = ex[a]*ey[a];

    //cc_xx[a] = ex[a]*ex[a];
    //cc_yy[a] = ey[a]*ey[a];
    //cc_xy[a] = ex[a]*ey[a];
	}

	for (ts=0; ts<=time; ts++)
	{
		for (i=1; i<=nx; i++)
		{
			for (j=0; j<=ny+1; j++)
			{
				/*tmp1 = pow((pow(i-xc,2.0) + pow(j-yc,2.0)),0.5);

				if (tmp1 <= 0.5*D)	isn[i][j] = 1;
				else								isn[i][j] = 0;*/
        isn[i][j] = 0;
			  if (j==0   ) {isn[i][j] = 1;}
			  if (j==ny+1) {isn[i][j] = 1;}
			}
		}
		//----------------------------------------------------------------------
		rhoAvg = 0.0;

		for (i=1; i<=nx; i++)
		{
			for (j=1; j<=ny; j++)
			{
				tmp1 = 0.0;
				tmp2 = 0.0;
				tmp3 = 0.0;

				for (a=0; a<q; a++)
				{
					tmp1 += f[a][i][j];
					tmp2 += f[a][i][j]*ex[a];
					tmp3 += f[a][i][j]*ey[a];
				}

				rho[i][j] = tmp1;
				ux[i][j]	= tmp2/tmp1;
				uy[i][j]	= tmp3/tmp1;

				rhoAvg	 += tmp1;
			}
		}

		rhoAvg /=(nx*ny);
		//----------------------------------------------------------------------
		for (i=1; i<=nx; i++)
		{
			for (j=1; j<=ny; j++)
			{
				for (a=0; a<q; a++)
				{
					cDOTu				= ux[i][j]*ex[a]  + uy[i][j]*ey[a] ;
					uDOTu				= ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
					feq			= wt[a]*rho[i][j]*(1.0 + 3.0*cDOTu + 4.5*cDOTu*cDOTu - 1.5*uDOTu);
					ft[a][i][j]	= f[a][i][j] - (f[a][i][j]-feq)/tau; //collision
				}
			}
		}
		//----------------------------------------------------------------------
		for (i=1; i<=nx; i++) //Streaming post-collision
		{
			for (j=1; j<=ny; j++)
			{
				for (a=0; a<q; a++)
				{
					ia = i+ex[a];
					ja = j+ey[a];

					//if (ia<1 )	{ ia = nx;  }
					//if (ia>nx)	{ ia = 1 ;  }

					f[a][ia][ja] = ft[a][i][j];
				}
			}
		}
		//----------------------------------------------------------------------
		/*for(j=1,i=1;i<=nx;i++)
		{
			f[5][i][j]=f[8][i][j-1];
			f[2][i][j]=f[4][i][j-1];
			f[6][i][j]=f[7][i][j-1];
		}

		for(j=ny,i=1;i<=nx;i++)
		{
			f[8][i][j]=f[5][i][j+1];
			f[4][i][j]=f[2][i][j+1];
			f[7][i][j]=f[6][i][j+1];
		}*/

for(j=1;j<=ny;j++)
		{
			//Inlet
			i=1;
			ux[i][j] = uInlet;//6*uInlet*(ny-(j-0.5))*(j-0.5)/(ny*ny);
			uy[i][j] = 0.0;
			rho[i][j] = (f[0][i][j] + f[2][i][j] + f[4][i][j] + 2.0*(f[3][i][j] + f[6][i][j] + f[7][i][j]))/(1-ux[i][j]);
      //------term1-------
			dudx = (-3.0*ux[i][j] + 4.0*ux[i+1][j] - 1.0*ux[i+2][j])/2.0; //dudx
			dvdx = (-3.0*uy[i][j] + 4.0*uy[i+1][j] - 1.0*uy[i+2][j])/2.0;
      if(j == 0)
			{
				dudy = (-3.0*ux[i][j] + 4.0*ux[i][j+1] - 1.0*ux[i][j+2])/2.0;
				dvdy = (-3.0*uy[i][j] + 4.0*uy[i][j+1] - 1.0*uy[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				dudy = (3.0*ux[i][j] - 4.0*ux[i][j-1] + 1.0*ux[i][j-2])/2.0;
				dvdy = (3.0*uy[i][j] - 4.0*uy[i][j-1] + 1.0*uy[i][j-2])/2.0;
			}
			else
			{
				dudy = (ux[i][j+1] - ux[i][j-1])/2.0;
				dvdy = (uy[i][j+1] - uy[i][j-1])/2.0;
			}
			//-----term2-----
			drhou2dx = (-3.0*rho[i][j]*ux[i][j]*ux[i][j] + 4.0*rho[i+1][j]*ux[i+1][j]*ux[i+1][j] - 1.0*rho[i+2][j]*ux[i+2][j]*ux[i+2][j])/2.0; //dudx
			drhouvdx = (-3.0*rho[i][j]*ux[i][j]*uy[i][j] + 4.0*rho[i+1][j]*ux[i+1][j]*uy[i+1][j] - 1.0*rho[i+2][j]*ux[i+2][j]*uy[i+2][j])/2.0; //dudx
			if(j == 0)
			{
				drhouvdy = (-3.0*rho[i][j]*ux[i][j]*uy[i][j] + 4.0*rho[i][j+1]*ux[i][j+1]*uy[i][j+1] - 1.0*rho[i][j+2]*ux[i][j+2]*uy[i][j+2])/2.0;
				drhovvdy = (-3.0*rho[i][j]*uy[i][j]*uy[i][j] + 4.0*rho[i][j+1]*uy[i][j+1]*uy[i][j+1] - 1.0*rho[i][j+2]*uy[i][j+2]*uy[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				drhouvdy = (-3.0*rho[i][j]*ux[i][j]*uy[i][j] + 4.0*rho[i][j-1]*ux[i][j-1]*uy[i][j-1] - 1.0*rho[i][j-2]*ux[i][j-2]*uy[i][j-2])/2.0;
				drhovvdy = (-3.0*rho[i][j]*uy[i][j]*uy[i][j] + 4.0*rho[i][j-1]*uy[i][j-1]*uy[i][j-1] - 1.0*rho[i][j-2]*uy[i][j-2]*uy[i][j-2])/2.0;
			}
			else
			{
				drhouvdy = (rho[i][j]*ux[i][j+1]*uy[i][j+1] - rho[i][j]*ux[i][j-1]*uy[i][j-1])/2.0;
				drhovvdy = (rho[i][j]*uy[i][j+1]*uy[i][j+1] - rho[i][j]*uy[i][j-1]*uy[i][j-1])/2.0;
			}

			uDOTu = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
			for(a=0;a<9;a++)
			{
				//------term3-------
				QddRhouu   =Q_xx[a]*rho[i  ][j]*ux[i  ][j]*ux[i  ][j] + 2.0*Q_xy[a]*rho[i  ][j]*ux[i  ][j]*uy[i  ][j] + Q_yy[a]*rho[i  ][j]*uy[i  ][j]*uy[i  ][j];
        QddRhouuxp1=Q_xx[a]*rho[i+1][j]*ux[i+1][j]*ux[i+1][j] + 2.0*Q_xy[a]*rho[i+1][j]*ux[i+1][j]*uy[i+1][j] + Q_yy[a]*rho[i+1][j]*uy[i+1][j]*uy[i+1][j];
        QddRhouuxp2=Q_xx[a]*rho[i+2][j]*ux[i+2][j]*ux[i+2][j] + 2.0*Q_xy[a]*rho[i+2][j]*ux[i+2][j]*uy[i+2][j] + Q_yy[a]*rho[i+2][j]*uy[i+2][j]*uy[i+2][j];
        dQddRhouudx=(-3.0*QddRhouu + 4.0*QddRhouuxp1 - 1.0*QddRhouuxp2)/2.0;

        QddRhouuyp1=Q_xx[a]*rho[i][j+1]*ux[i][j+1]*ux[i][j+1] + 2.0*Q_xy[a]*rho[i][j+1]*ux[i][j+1]*uy[i][j+1] + Q_yy[a]*rho[i][j+1]*uy[i][j+1]*uy[i][j+1];
        QddRhouuym1=Q_xx[a]*rho[i][j-1]*ux[i][j-1]*ux[i][j-1] + 2.0*Q_xy[a]*rho[i][j-1]*ux[i][j-1]*uy[i][j-1] + Q_yy[a]*rho[i][j-1]*uy[i][j-1]*uy[i][j-1];
   			if(j == 0)
  			{
          QddRhouuyp2=Q_xx[a]*rho[i][j+2]*ux[i][j+2]*ux[i][j+2] + 2.0*Q_xy[a]*rho[i][j+2]*ux[i][j+2]*uy[i][j+2] + Q_yy[a]*rho[i][j+2]*uy[i][j+2]*uy[i][j+2];
  				dQddRhouudy=(-3.0*QddRhouu + 4.0*QddRhouuyp1 - 1.0*QddRhouuyp2)/2.0;
  			}
  			else if(j == ny)
  			{
          QddRhouuym2=Q_xx[a]*rho[i][j-2]*ux[i][j-2]*ux[i][j-2] + 2.0*Q_xy[a]*rho[i][j-2]*ux[i][j-2]*uy[i][j-2] + Q_yy[a]*rho[i][j-2]*uy[i][j-2]*uy[i][j-2];
  				dQddRhouudy = (3.0*QddRhouu - 4.0*QddRhouuym1 + 1.0*QddRhouuym2)/2.0;
  			}
  			else
  			{
  				dQddRhouudy = (QddRhouuyp1 - QddRhouuym1)/2.0;
  			}
				term1 = -wt[a]*cs2Inv*tau*rho[i][j]*(Q_xx[a]*dudx + Q_xy[a]*dudy + Q_xy[a]*dvdx + Q_yy[a]*dvdy);
				term2 = -wt[a]*cs2Inv*tau*(ex[a]*drhou2dx + ex[a]*drhouvdy + ey[a]*drhouvdx +ey[a]*drhovvdy);
				term3 = -wt[a]*cs2Inv*cs2Inv*tau*0.5*(ex[a]*dQddRhouudx + ey[a]*dQddRhouudy);
				f_neq = term1 - term2 + term3;
				cDOTu = ex[a]*ux[i][j] + ey[a]*uy[i][j];
				feq = wt[a]*rho[i][j]*(1 + cDOTu*cs2Inv + 0.5*cDOTu*cDOTu*cs2Inv*cs2Inv - 0.5*cs2Inv*uDOTu);
				f[a][i][j] = feq + f_neq;
			}

			//Outlet
			i = nx;
			rho[i][j] = 1.0;
		 	ux[i][j]=(f[0][i][j]+f[2][i][j]+f[4][i][j]+2*(f[1][i][j]+f[5][i][j]+f[8][i][j]))/rho[i][j] - 1;
      uy[i][j]=0.0;
			//------term1-------
			dudx = -(-3.0*ux[i][j] + 4.0*ux[i-1][j] - 1.0*ux[i-2][j])/2.0; //dudx
			dvdx = -(-3.0*uy[i][j] + 4.0*uy[i-1][j] - 1.0*uy[i-2][j])/2.0;
			if(j == 0)
			{
				dudy = (-3.0*ux[i][j] + 4.0*ux[i][j+1] - 1.0*ux[i][j+2])/2.0;
				dvdy = (-3.0*uy[i][j] + 4.0*uy[i][j+1] - 1.0*uy[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				dudy = (3.0*ux[i][j] - 4.0*ux[i][j-1] + 1.0*ux[i][j-2])/2.0;
				dvdy = (3.0*uy[i][j] - 4.0*uy[i][j-1] + 1.0*uy[i][j-2])/2.0;
			}
			else
			{
				dudy = (ux[i][j+1] - ux[i][j-1])/2.0;
				dvdy = (uy[i][j+1] - uy[i][j-1])/2.0;
			}
			//-----term2-----
			drhou2dx = -(-3.0*rho[i][j]*ux[i][j]*ux[i][j] + 4.0*rho[i-1][j]*ux[i-1][j]*ux[i-1][j] - 1.0*rho[i-2][j]*ux[i-2][j]*ux[i-2][j])/2.0; //dudx
			drhouvdx = -(-3.0*rho[i][j]*ux[i][j]*uy[i][j] + 4.0*rho[i-1][j]*ux[i-1][j]*uy[i-1][j] - 1.0*rho[i-2][j]*ux[i-2][j]*uy[i-2][j])/2.0; //dudx
			if(j == 0)
			{
				drhouvdy = (-3.0*rho[i][j]*ux[i][j]*uy[i][j] + 4.0*rho[i][j+1]*ux[i][j+1]*uy[i][j+1] - 1.0*rho[i][j+2]*ux[i][j+2]*uy[i][j+2])/2.0;
				drhovvdy = (-3.0*rho[i][j]*uy[i][j]*uy[i][j] + 4.0*rho[i][j+1]*uy[i][j+1]*uy[i][j+1] - 1.0*rho[i][j+2]*uy[i][j+2]*uy[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				drhouvdy = (-3.0*rho[i][j]*ux[i][j]*uy[i][j] + 4.0*rho[i][j-1]*ux[i][j-1]*uy[i][j-1] - 1.0*rho[i][j-2]*ux[i][j-2]*uy[i][j-2])/2.0;
				drhovvdy = (-3.0*rho[i][j]*uy[i][j]*uy[i][j] + 4.0*rho[i][j-1]*uy[i][j-1]*uy[i][j-1] - 1.0*rho[i][j-2]*uy[i][j-2]*uy[i][j-2])/2.0;
			}
			else
			{
				drhouvdy = (rho[i][j]*ux[i][j+1]*uy[i][j+1] - rho[i][j]*ux[i][j-1]*uy[i][j-1])/2.0;
				drhovvdy = (rho[i][j]*uy[i][j+1]*uy[i][j+1] - rho[i][j]*uy[i][j-1]*uy[i][j-1])/2.0;
			}

			uDOTu = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
			for(a=0;a<9;a++)
			{
				//------term3-------
				QddRhouu   =Q_xx[a]*rho[i  ][j]*ux[i  ][j]*ux[i  ][j] + 2.0*Q_xy[a]*rho[i  ][j]*ux[i  ][j]*uy[i  ][j] + Q_yy[a]*rho[i  ][j]*uy[i  ][j]*uy[i  ][j];
				QddRhouuxp1=Q_xx[a]*rho[i-1][j]*ux[i-1][j]*ux[i-1][j] + 2.0*Q_xy[a]*rho[i-1][j]*ux[i-1][j]*uy[i-1][j] + Q_yy[a]*rho[i-1][j]*uy[i-1][j]*uy[i-1][j];
				QddRhouuxp2=Q_xx[a]*rho[i-2][j]*ux[i-2][j]*ux[i-2][j] + 2.0*Q_xy[a]*rho[i-2][j]*ux[i-2][j]*uy[i-2][j] + Q_yy[a]*rho[i-2][j]*uy[i-2][j]*uy[i-2][j];
				dQddRhouudx=-(-3.0*QddRhouu + 4.0*QddRhouuxp1 - 1.0*QddRhouuxp2)/2.0;

				QddRhouuyp1=Q_xx[a]*rho[i][j+1]*ux[i][j+1]*ux[i][j+1] + 2.0*Q_xy[a]*rho[i][j+1]*ux[i][j+1]*uy[i][j+1] + Q_yy[a]*rho[i][j+1]*uy[i][j+1]*uy[i][j+1];
				QddRhouuym1=Q_xx[a]*rho[i][j-1]*ux[i][j-1]*ux[i][j-1] + 2.0*Q_xy[a]*rho[i][j-1]*ux[i][j-1]*uy[i][j-1] + Q_yy[a]*rho[i][j-1]*uy[i][j-1]*uy[i][j-1];
				if(j == 0)
				{
					QddRhouuyp2=Q_xx[a]*rho[i][j+2]*ux[i][j+2]*ux[i][j+2] + 2.0*Q_xy[a]*rho[i][j+2]*ux[i][j+2]*uy[i][j+2] + Q_yy[a]*rho[i][j+2]*uy[i][j+2]*uy[i][j+2];
					dQddRhouudy=(-3.0*QddRhouu + 4.0*QddRhouuyp1 - 1.0*QddRhouuyp2)/2.0;
				}
				else if(j == ny)
				{
					QddRhouuym2=Q_xx[a]*rho[i][j-2]*ux[i][j-2]*ux[i][j-2] + 2.0*Q_xy[a]*rho[i][j-2]*ux[i][j-2]*uy[i][j-2] + Q_yy[a]*rho[i][j-2]*uy[i][j-2]*uy[i][j-2];
					dQddRhouudy = (3.0*QddRhouu - 4.0*QddRhouuym1 + 1.0*QddRhouuym2)/2.0;
				}
				else
				{
					dQddRhouudy = (QddRhouuyp1 - QddRhouuym1)/2.0;
				}
				term1 = -wt[a]*cs2Inv*tau*rho[i][j]*(Q_xx[a]*dudx + Q_xy[a]*dudy + Q_xy[a]*dvdx + Q_yy[a]*dvdy);
				term2 = -wt[a]*cs2Inv*tau*(ex[a]*drhou2dx + ex[a]*drhouvdy + ey[a]*drhouvdx +ey[a]*drhovvdy);
				term3 = -wt[a]*cs2Inv*cs2Inv*tau*0.5*(ex[a]*dQddRhouudx + ey[a]*dQddRhouudy);
				f_neq = term1 - term2 + term3;
				cDOTu = ex[a]*ux[i][j] + ey[a]*uy[i][j];
				feq = wt[a]*rho[i][j]*(1 + cDOTu*cs2Inv + 0.5*cDOTu*cDOTu*cs2Inv*cs2Inv - 0.5*cs2Inv*uDOTu);
				f[a][i][j] = feq + f_neq;
			}
		}
		 /*for(i=1,j=1;j<=ny;j++)
		 {
		 	rho[i][j]=(f[0][i][j]+f[2][i][j]+f[4][i][j]+2*(f[6][i][j]+f[3][i][j]+f[7][i][j]))/(1-uInlet);
		 	f[1][i][j]=f[3][i][j]+((2.0/3.0)*rho[i][j]*uInlet);
		 	f[5][i][j]=f[7][i][j]-(0.5*(f[2][i][j]-f[4][i][j]))+((1.0/6.0)*rho[i][j]*uInlet);
		 	f[8][i][j]=f[6][i][j]+(0.5*(f[2][i][j]-f[4][i][j]))+((1.0/6.0)*rho[i][j]*uInlet);
		 }

		 for(i=nx,j=1;j<=ny;j++)
		 {
		 	ux[i][j]=(f[0][i][j]+f[2][i][j]+f[4][i][j]+2*(f[1][i][j]+f[5][i][j]+f[8][i][j]))/rho0 - 1;
		 	f[3][i][j]=f[1][i][j]-((2.0/3.0)*rho0*ux[i][j]);
		 	f[6][i][j]=f[8][i][j]-0.5*(f[2][i][j]-f[4][i][j])-((1.0/6.0)*rho0*ux[i][j]);
		 	f[7][i][j]=f[5][i][j]+0.5*(f[2][i][j]-f[4][i][j])-((1.0/6.0)*rho0*ux[i][j]);
		 }*/
		//----------------------------------------------------------------------

		for (i=1; i<=nx; i++) //BC
		{
			for (j=1; j<=ny; j++)
			{
				if (isn[i][j] == 0)
				{
					for (a=0; a<q; a++)
					{
						ia = i+ex[a];
						ja = j+ey[a];

            if (isn[ia][ja] == 1) //Wall
						{
							f[kb[a]][i][j] = f[a][ia][ja];
						}

						/*if (isn[ia][ja]==1) //Cylinder
						{
							f[kb[a]][i ][j ] = ft[a		 ][i ][j ];
							f[a    ][ia][ja] = ft[kb[a]][ia][ja];

							fx[1] += ex[a]*2.0*(-ft[kb[a]][ia][ja] + ft[a][i][j]);
							fy[1] += ey[a]*2.0*(-ft[kb[a]][ia][ja] + ft[a][i][j]);
						}*/
					}
				}
			}
		}
		//----------------------------------------------------------------------
		//fx_t 	= 0.5*(fx[0]  +  fx[1]);
		//fy_t 	= 0.5*(fy[0]  +  fy[1]);

		//Cd = fx_t/(0.5*rho0*uInlet*uInlet*D);
		//Cl = fy_t/(0.5*rho0*uInlet*uInlet*D);
		//----------------------------------------------------------------------
		if (ts % dispFreq == 0)
		{
			fprintf(fid,"%d\t%10.6f\t%10.6f\n", ts, rhoAvg, ux[9*nx/10][(ny+1)/2]/(1.5*uInlet));
			printf(     "%d\t%10.6f\n",	ts, rhoAvg);
		}
		//----------------------------------------------------------------------
		//fx[0] = fx[1];
		//fy[0] = fy[1];
		//----------------------------------------------------------------------
		if(ts<=time && ts % (time/(noOfSnaps-1)) == 0)
		{
			solnumber++;

			strcpy(filename,prefix);
			sprintf(solstr,"%d",solnumber);
			strcat(filename,solstr);
			strcat(filename,type);
			soln = fopen(filename,"w");

			fprintf(soln,"Variables=x,y,u,v,rho,region\n");
			fprintf(soln,"Zone I= %d,J= %d\n\n",nx,ny);

			for (j=1; j<=ny; j++)
			{
				for (i=1; i<=nx; i++)
				{
					fprintf(soln,"%d %d %12.8f %12.8f %12.8f %d\n",i,j,ux[i][j],uy[i][j],rho[i][j],isn[i][j]);
				}
				fprintf(soln,"\n");
			}
			fclose(soln);
			printf("snap %d recorded at time %d\n",solnumber,ts);
		}
		//----------------------------------------------------------------------
	}//Time loop Ends

  for (i=1; i<=nx; i++)
	{
    fprintf(fid2,"%d\t%10.6f\t%10.6f\t%10.6f\n", i, rho[i][(ny+1)/2], ux[i][(ny+1)/2], uy[i][(ny+1)/2]);
   }

   for (j=1; j<=ny; j++)
	{
    fprintf(fid3,"%d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n", j, ux_a[j]/(1.5*uInlet), ux[nx/4][j]/(1.5*uInlet), ux[nx/2][j]/(1.5*uInlet), ux[9*nx/10][j]/(1.5*uInlet));
   }

	fclose(fid);
  fclose(fid2);
  fclose(fid3);
	return 0;
}
