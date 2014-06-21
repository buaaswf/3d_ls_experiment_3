#include "ThreeDim_LevelSet.h"
void peppersalt(Raw *data,double percent);
//one ball
Raw * myColondata1(int l,int m, int n)
{
	int R=100;
	Raw * data= new Raw(l,m,n);

	//memset(data->getdata(),0,l*m*n);

	
	for (int i = 0 ; i< l; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < n; k++)
			{
				if ((i-l/2)*(i-l/2) + (j-m/2)*(j-m/2) + (k-n/2)*(k-n/2)  < R*R)
				{
					data->put(i,j,k,1000);
				}
				else if ((i-l/2)*(i-l/2) + (j-m/2)*(j-m/2) + (k-n/2)*(k-n/2) < (R+10)*(R+10))
				{
					data->put(i,j,k,0);
				}else if (data->get(i,j,k)!=2000)
				{
					data->put(i,j,k,-1000);
				}





			}
		}
	}
	peppersalt(data,0.99);
	return data;
}
//two balls
Raw * myColondata2(int l,int m, int n)
{
	int R=115;
	Raw * data= new Raw(l,m,n);

	//memset(data->getdata(),0,l*m*n);


	for (int i = 0 ; i< l; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < n; k++)
			{
				if ((i-l/4)*(i-l/4) + (j-m*2/4)*(j-m*2/4) + (k-n/2)*(k-n/2)  < R*R||
					(i-l*3/4)*(i-l*3/4) + (j-m*2/4)*(j-m*2/4) + (k-n/2)*(k-n/2)  < R*R)
				{
					data->put(i,j,k,1000);
				}
				else if ((i-l/4)*(i-l/4) + (j-m*2/4)*(j-m*2/4) + (k-n/2)*(k-n/2)  < (R+10)*(R+10)
					||(i-l*3/4)*(i-l*3/4) + (j-m*2/4)*(j-m*2/4) + (k-n/2)*(k-n/2) < (R+10)*(R+10))
				{
					data->put(i,j,k,0);
				}else if (data->get(i,j,k)!=2000)
				{
					data->put(i,j,k,-1000);
				}





			}
		}
	}
	peppersalt(data,0.99);
	return data;
}
Raw * myColondata3(int l,int m, int n)
{
	int R=50;
	int yR=50;
	Raw * data= new Raw(l,m,n);

	//memset(data->getdata(),0,l*m*n);
	int zR=115;
	int maxr=10;
	for (int i = 0 ; i< l; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < n; k++)
			{
				for (int r=0; r<maxr; r++)
				{
					//if (yR*yR*zR*zR*(i-l/4)*(i-l/4) +zR*zR*R*R*(j-m*r/maxr)*(j-m*r/maxr) + yR*yR*R*R*(k-n/2)*(k-n/2)  < zR*zR*yR*yR*R*R||
					//	yR*yR*zR*zR*(i-l*3/4)*(i-l*3/4) + zR*zR*R*R*(j-m*r/maxr)*(j-m*r/maxr) + yR*yR*R*R*(k-n/2)*(k-n/2)  < zR*zR*yR*yR*R*R)
					//{
					//	data->put(i,j,k,1000);
					//}
					// 
					if ((i-l/4)*(i-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R||
						(i-l*2/4)*(i-l*2/4) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R)
					{
						data->put(i,j,k,1000);
					}
					else if ((i-l/4)*(i-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R+10)*(R+10)||
						(i-l*2/4)*(i-l*2/4) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R+10)*(R+10))
					{
						if (data->get(i,j,k)!=1000)
						{
							data->put(i,j,k,0);
						}
						
					}
					else if (data->get(i,j,k)!=2000&& data->get(i,j,k)!=1000&&data->get(i,j,k)!=0)
					{
						data->put(i,j,k,-1000);
					}
				}
		





			}
		}
	}
	peppersalt(data,0.99);
	return data;
}
Raw * myColondata4(int l,int m, int n,double percent)
{
	int R=50;
	int yR=50;
	Raw * data= new Raw(l,m,n);

	//memset(data->getdata(),0,l*m*n);
	int zR=115;
	int maxr=10;
	int outerthickness= 15;
	int distance=10;
	for (int i = 0 ; i< l; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < n; k++)
			{
				for (int r=0; r<maxr; r++)
				{
					if (r<5)
					{
						if ((i-sqrt((double)m*r/maxr)-l/4)*(i-sqrt((double)m*r/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R)			//the || represents two colons
						{
							data->put(i,j,k,1000);
						}
						else if ((i-sqrt((double)m*r/maxr)-l/4)*(i-sqrt((double)m*r/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R + outerthickness)*(R + outerthickness)||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R + outerthickness)*(R + outerthickness))   //the || represents two colons
						{
							if (data->get(i,j,k)!=1000)
							{
								data->put(i,j,k,0);
							}

						}
						else if (data->get(i,j,k)!=2000&& data->get(i,j,k)!=1000&&data->get(i,j,k)!=0)
						{
							data->put(i,j,k,-1000);
						}

					}
					else 
					{
						//r-=5;
						if ((i+sqrt((double)m*(r-5)/maxr)-l/4)*(i+sqrt((double)m*(r-5)/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R)
						{
							data->put(i,j,k,1000);
						}
						else if ((i+sqrt((double)m*(r-5)/maxr)-l/4)*(i+sqrt((double)m*(r-5)/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R + outerthickness)*(R + outerthickness)||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R+ outerthickness)*(R + outerthickness))
						{
							if (data->get(i,j,k)!=1000)
							{
								data->put(i,j,k,0);
							}

						}
						else if (data->get(i,j,k)!=2000&& data->get(i,j,k)!=1000&&data->get(i,j,k)!=0)
						{
							data->put(i,j,k,-1000);
						}
						//r+=5;
					}
				}






			}
		}
	}
	peppersalt(data,percent);
	return data;
}
Raw * myColondata5(int l,int m, int n,double percent)
{
	int R=50;
	int yR=50;
	Raw * data= new Raw(l,m,n);

	//memset(data->getdata(),0,l*m*n);
	int zR=115;
	int maxr=10;
	int outerthickness= 15;
	int distance=10;
	for (int i = 0 ; i< l; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < n; k++)
			{
				for (int r=0; r<maxr; r++)
				{
					if (r<5)
					{
						if ((i-sqrt((double)m*r/maxr)-l/4)*(i-sqrt((double)m*r/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R)			//the || represents two colons
						{
							data->put(i,j,k,1000);
						}
						else if ((i-sqrt((double)m*r/maxr)-l/4)*(i-sqrt((double)m*r/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R + outerthickness)*(R + outerthickness)||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R + outerthickness)*(R + outerthickness))   //the || represents two colons
						{
							if (data->get(i,j,k)!=1000)
							{
								data->put(i,j,k,-(i-sqrt((double)m*r/maxr)-l/4)*(i-sqrt((double)m*r/maxr)-l/4));
							}

						}
						else if (data->get(i,j,k)!=2000&& data->get(i,j,k)!=1000&&data->get(i,j,k)>0)
						{
							data->put(i,j,k,-1000);
						}

					}
					else 
					{
						//r-=5;
						if ((i+sqrt((double)m*(r-5)/maxr)-l/4)*(i+sqrt((double)m*(r-5)/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < R*R)
						{
							data->put(i,j,k,1000);
						}
						else if ((i+sqrt((double)m*(r-5)/maxr)-l/4)*(i+sqrt((double)m*(r-5)/maxr)-l/4) +(j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R + outerthickness)*(R + outerthickness)||
							(i-l*2/4-distance)*(i-l*2/4-distance) + (j-m*r/maxr)*(j-m*r/maxr) + (k-n/2)*(k-n/2)  < (R+ outerthickness)*(R + outerthickness))
						{
							if (data->get(i,j,k)!=1000)
							{
								data->put(i,j,k,-(i-sqrt((double)m*r/maxr)-l/4)*(i-sqrt((double)m*r/maxr)-l/4));
							}

						}
						else if (data->get(i,j,k)!=2000&& data->get(i,j,k)!=1000&&data->get(i,j,k)>0)
						{
							data->put(i,j,k,-1000);
						}
						//r+=5;
					}
				}






			}
		}
	}
	peppersalt(data,percent);
	return data;
}
void peppersalt(Raw *data, double percent)
{
	float pepper=0.5;
	float salt=0.5;
	double temp;
	double SNR=percent;
	    /*-----------------------------------------------------------------------------
     *  Salt & Pepper should be between 0 and 1.
     *-----------------------------------------------------------------------------*/
    //salt -= int (salt);
    //pepper -= int (pepper);
    /*-----------------------------------------------------------------------------
     *  Init the random seed with time.
     *-----------------------------------------------------------------------------*/
    srand (time (NULL));
    /*-----------------------------------------------------------------------------
     *  Salt & Pepper.
     *-----------------------------------------------------------------------------*/
	int size= (int)(data->size() * (1-SNR));
	srand(time(NULL));
	for(int i=0; i<size; i++) {

		int row = (int)(rand()%data->getXsize());
		int col = (int)(rand()%data->getYsize());
		int zsize=(int)(rand()%data->getZsize());
		data->put(row,col,zsize,2000);

	}
}