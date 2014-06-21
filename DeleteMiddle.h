
/************************************************************************/
/* need to impl the distance transform                                  */
/*********************************************************************/
//the levleset thickness data and the origion data threshold .the levleset data is smooth,
// but the middle of the two different regions may be growing into one region,which is not the truth.
// the threshold method can grow into two different regions.but the geometry is not good and not either smooth.
// so put the two method together.
//first,the origion data.raw l,m,n;
// second the thickness data after the  distance transform
// condition:if the thickness data's value is larger than the threshold(normally,3),
// && the position's value of the origin data is in the[-50,+wuqiongda],then put the thickness data the position zero.
// output the new thickness data.

void deleteMiddle(int argc,string dir)
{

	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=10;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\"
	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	strcat(dirhead,dirbody);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream("K:\\sdf\\volume\\clean\\clean\\ep\\clean\\m\\WI_3035_P_iso_clean.raw",&l,&m,&n);
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	for (int i = 0; i < l*m*n; i++)
	{
		inputo[i]=(float) indata[i];		
	}

	Raw *input=new Raw(l,m,n,inputo);
	delete [] indata;
	float *thicnknessdistance=test.readStreamfloat("D:\\swfdata20140420res\\58_3035P\\edt\\3035p-littleedt.raw",&l,&m,&n);//read the data after distance transform
	unsigned char *thicnknessorigin = new unsigned char[l*m*n];
	test.readImage(thicnknessorigin,"D:\\swfdata20140420res\\58_3035P\\thickness\\WI_3035_P_iso_clean.rawthickness5-8_2_20140405.raw",l*m*n);
	float *thicknessf= new float[l*m*n];
	for (int i=0;i<l*m*n;i++)
	{
		thicknessf[i]=(float(thicnknessorigin[i]));
	}
	Raw * thickness = new Raw(l,m,n,thicnknessdistance);
	Raw *originthickness = new Raw(l,m,n,thicknessf);
	
	PIXTYPE maxTissue= -100;
	for (int i = 0; i < thickness->getXsize(); i++)
	{
		for (int j=0; j <thickness->getYsize(); j++)
		{
			for ( int k=0; k < thickness->getZsize(); k++)
			{
				if (thickness->get(i,j,k)>0)
				{
				
						int halfsize=2;
						int countneighbour=0;
						PIXTYPE val = input->get(i,j,k); 
						//PIXTYPE val2 = originthickness->get(i,j,k); 
						//for (int l=i-halfsize;l<i+halfsize;l++)
						//{
						//	for (int m=j-halfsize; m<j+halfsize;m++)
						//	{
						//		for (int n=k-halfsize; n<k+halfsize;n++)
						//		{
						//			if (val > maxTissue)
						//			{
						//				countneighbour++;
						//			}
						//		}
						//	}
						//}
						//if (countneighbour>10 )
						//{
						//	originthickness->put(k,j,i,0);
						//}
						if ( val>= -124 || val <= -924)
						{
							originthickness->put(i,j,k,0);
						}

					}//...if thickness >2
			
			}
		}
	}
	RawImage *write=new RawImage();
	write->writeImagecolon(*originthickness);


}
