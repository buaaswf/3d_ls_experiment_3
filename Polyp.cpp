#include <iostream>
#include <fstream>
#include <string>
#include "Polyp.h"
#include <queue>
using namespace std;
void seedGrowing(Raw * src,vector<Seed> seedlist);
void seedGrowing(Raw * src, Raw *origin, vector<Seed> seedlist);
void seedGrowing(Raw * src, vector<Seed> seedlist, int threshold);
Polyp::Polyp(void)
{
}


Polyp::~Polyp(void)
{
}
int CountLines(char *filename)
{
	ifstream ReadFile;
	int n=0;
	string tmp;
	ReadFile.open(filename,ios::in);//ios::in :readonly
	if(ReadFile.fail())// open fails
	{
		return 0;
	}
	else//if exit
	{
		while(getline(ReadFile,tmp,'\n'))
		{
			n++;
		}
		ReadFile.close();
		return n;
	}
}
vector<Seed> readSeedFromTXT(char *txtfilepath)
{
		ifstream file;
		int LINES;
		char filename[512]="inFile.txt";
		file.open(filename,ios::in);
		if(file.fail())
		{
			cout<<"not exit."<<endl;
			file.close();
		}
		else//if exit
		{
			LINES=CountLines(filename);
			int *tempInt=new int[LINES];
			char *tempChar=new char[LINES];
			int i=0;
			while(!file.eof()) //read data to array
			{

				file>>tempInt[i];
				file>>tempChar[i];
				i++;
			}
			file.close(); //close file
			for(i=0;i<LINES;i++)
				cout<<tempInt[i]<<"\t"<<tempChar[i]<<endl;
			delete []tempInt;
			delete []tempChar;
		}

		vector<Seed> seed;
		return seed;

}
Raw* Polyp::initialRegion(vector<Seed> seedlist,int size,int l, int m, int n)
{
	Raw *initialRegion=new Raw(l,m,n);
	//memset(initialRegion->getdata(),2,initialRegion->size());
	for (int i = 0; i < l*m*n; i++)
	{
		initialRegion->putXYZ(i,-2);
	}
	for (vector<Seed>::iterator it = seedlist.begin();it!=seedlist.end();++it)
	{
		
		//Seed *seed=new Seed(it->x,it->y,it->z);
		for (int i = it->x - size; i < it->x + size; ++i )
		{
			for (int j = it->y - size; j < it->y + size; ++j)
			{
				for (int k = it->z - size; k< it->z + size; ++k)
				{
					initialRegion->put(i,j,k,2);
				}
				
			}
			
		}
		
	}
	return initialRegion;
}
void Polyp::polypDetect(string dir)
{
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=20;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\

	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	strcat(dirhead,dirbody);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream(dirhead,&l,&m,&n);
	unsigned char * thicknessdata=test.readStreamuchar(dirhead,l,m,n,212);
	n=50;
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	for (int i = 0; i < l*m*n; i++)
	{
		inputo[i]=(float) indata[i];		
	}

	Raw *input=new Raw(l,m,n,inputo);
	/*Filter *f=new Filter();*/
	//input=f->guass3DFilter(input,3);
	RawImage *write=new RawImage();
	ThreeDim_LevelSet *ls=new ThreeDim_LevelSet();
	//20140405 delete because of the existance of 
	ls->initialg(*input);
	//for (int i=0; i<input->getXsize(); i++)
	//{
	//	for (int j=0; j<input->getYsize(); j++)
	//	{
	//		for (int k=0; k<input->getZsize(); k++)
	//		{
	//			if (input->get(i,j,k)>=1)
	//			{
	//				initial->put(i,j,k,-2);
	//			}
	//			else initial->put(i,j,k,2);

	//		}
	//	}

	//}
	vector<Seed> seedlist;
	Seed *seed=new Seed(164, 373, 20);
	seedlist.push_back(*seed);
	Raw *initialdata=initialRegion(seedlist,20,l,m,n);
	*initial=ls->minimal_surface(*initialdata,*input,5.0,0.1,3,1.5,1,iter_outer,pt);
	char *outname1="inner5-8_2.raw";
	char outdir[200]=output;

	strcat(outdir,dirbody);
	strcat(outdir,outname1);
	//test.readImage2(initial->getdata(),outdir,l*m*n);
	test.writeImageName(*initial,outdir);
	//Raw temp(*initial);
	//ls->outerwallauto(*initial,*input,10,0.1,-6,1.5,1,10,pt);
	//*initial -=temp;
	//char *outname2="outer5-8_2_20140405.raw";
	//char outdir2[200]=output;
	//strcat(outdir2,dirbody);
	//strcat(outdir2,outname2);
	//test.writeImageName(*initial,outdir2);
	//evaluate(dir,l,m,n);
}
vector<Seed> probabiltyouterpointleastcost(Raw *data, Seed *center,int range)
{
	int x = center->x;
	int y = center->y;
	int z = center->z;
	vector<Seed> datalist;
	for (int i = 0; i < data->getXsize();i++)
	{
		for (int j = 0; j < data->getYsize(); j++)
		{
			for (int k = 0; k < data->getZsize(); k++)
			{
				if ((i - x)*(i - x) + (j - y)*(j - y) + (k - z)*(k - z) <= range*range)
				{
					Seed *seed = new Seed(i, j, k);
					datalist.push_back(*seed);
				}
			}
		}
	}
	return datalist;


}
Seed * findLeastCostSeed(vector<Seed> list,Seed *center)
{
	int distance = 100000;
	int pos = 0;
	for (int i = 0; i < list.size(); i++)
	{
		int val = (list[i].x - center->x)*(list[i].x - center->x) + (list[i].y - center->y)*(list[i].y - center->y)
			+ (list[i].z - center->z)*(list[i].z - center->z);
		val < distance*distance ? distance = sqrt(val) : distance = distance;
		pos = i;
	}
	return &(list[pos]);
}
void Polyp::polypDetect(string dir,string dirthickness,int methodo)
{
	int offset = 220;
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=50;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\
	
	char dirheadt[200]=inputt; 
	char dirbodyt[100];
	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	strcpy(dirbodyt,dirthickness.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	cout <<"dirbody thickness"<<dirbodyt<<endl;
	strcat(dirhead,dirbody);
	strcat(dirheadt,dirbodyt);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream(dirhead,&l,&m,&n);

	
	unsigned char * thicknessdata=test.readStreamuchar(dirheadt,l,m,n,offset);//3036P 212 3033P 90 3033P peduncated:90 sessile1:70 sessile2:248,204,58:60
	//PIXTYPE * 
	n=50;
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];

	//for (int i=0; i<input->getXsize(); i++)
	//{
	//	for (int j=0; j<input->getYsize(); j++)
	//	{
	//		for (int k=0; k<input->getZsize(); k++)
	//		{
	//			if (input->get(i,j,k)>=1)
	//			{
	//				initial->put(i,j,k,-2);
	//			}
	//			else initial->put(i,j,k,2);

	//		}
	//	}

	//}
	vector<Seed> seedlist;
	Seed *seed=new Seed(248, 208, 15);//3036P :145, 379, 25 //248, 208, 15 3033P:peducated:247, 208, 15 sessile1: 246,222,4+63sessile2:248, 204, 18:
	seedlist.push_back(*seed);
	int method=methodo;
	if (method==1)
	{
		for (int i = 0; i < l*m*n; i++)
		{
			inputo[i] = (float)indata[i];
		}
		delete[] thicknessdata;
		Raw *input = new Raw(l, m, n, inputo);
		/*Filter *f=new Filter();*/
		//input=f->guass3DFilter(input,3);
		RawImage *write = new RawImage();
		ThreeDim_LevelSet *ls = new ThreeDim_LevelSet();
		//20140405 delete because of the existance of 
		ls->initialg(*input);
	Raw *initialdata=initialRegion(seedlist,2,l,m,n);
	*initial=ls->minimal_surface(*initialdata,*input,2*5.0,0.1,-3,1.5,1,iter_outer,pt);
	char *outname1="inner5-8_2.raw";
	char outdir[200]=output;

	strcat(outdir,dirbody);
	strcat(outdir,outname1);
	//test.readImage2(initial->getdata(),outdir,l*m*n);
	test.writeImageName(*initial,outdir);
	//delete[] input;
	//Raw temp(*initial);
	//ls->outerwallauto(*initial,*input,10,0.1,-6,1.5,1,10,pt);
	//*initial -=temp;
	//char *outname2="outer5-8_2_20140405.raw";
	//char outdir2[200]=output;
	//strcat(outdir2,dirbody);
	//strcat(outdir2,outname2);
	//test.writeImageName(*initial,outdir2);
	//evaluate(dir,l,m,n);
	}
	else if (method==2)
	{
		for (int i = 0; i < l*m*n; i++)
		{
			inputo[i] = (float)thicknessdata[i];
		}
		PIXTYPE *origin = new PIXTYPE[l*m*(n+offset)];
		for (int j = 0; j < l*m*(n+offset); j++)
		{
			origin[j] = (float)indata[j];
		}

		Raw *origindata = new Raw(l, m, offset+n, origin);
		//delete[] thicknessdata;
		Raw *input = new Raw(l, m, n, inputo);
		/*Filter *f=new Filter();*/
		//input=f->guass3DFilter(input,3);
		RawImage *write = new RawImage();
		ThreeDim_LevelSet *ls = new ThreeDim_LevelSet();
		//20140405 delete because of the existance of 
		ls->initialg(*input);
		seedGrowing(input,seedlist);//fisrt time delete the 
		for (int i = 0; i<input->getXsize(); i++)
		{
			for (int j = 0; j<input->getYsize(); j++)
			{
				for (int k = 0; k<input->getZsize(); k++)
				{
					PIXTYPE val = input->get(i, j, k);
					//if (((k-256)*(k-256)+(j-256)*(j-256) )<(230*230))//k<409 && k> 107 && j>156 &&j <390
					//{
					if (val !=1)
					{
						val = 200;  //change to 100 for roc computing *val=0; 
						

					}
					else val = 0; ////change to 0 for roc computing *val=100; 
					input->put(i, j, k, val);
					//}
					//else *val = 0;
				}
			}
		}
		Seed *seed3 = new Seed(169,367,11);//250,208,12
		seedlist.clear();
		//seedlist.push_back(*seed3);
		//seedGrowing(input,seedlist,0);
		//Seed *seed = new Seed(248, 209, 18);//3036P :145, 379, 25 //248, 208, 15 3033P:peducated:247, 208, 15 sessile1: 246,222,4+63sessile2:248, 204, 18:
		//seedlist.push_back(*seed);
		//seedGrowing(input,seedlist,0);
		//ls->initialg(*input);
		seedlist.clear();
		Seed *seed2 = new Seed(169, 367, 11);//3036P :145, 379, 25 //248, 208, 15 3033P:peducated:247, 208, 15 sessile1: 246,222,4+63sessile2:248, 204, 18://250, 208, 12
		seedlist.push_back(*seed2);
		Raw *initialdata = initialRegion(seedlist,5, l, m, n);
		*initial = ls->minimal_surface(*initialdata, *input, 1.0, 0.1, 1, 1.5, 1, iter_outer, pt);//alfa is the ballon pressure,lambda control the curvature
		char *outname1="inner5-8_2polypmethod2.raw";
		char outdir[200]=output;
		strcat(outdir,dirbody);
		strcat(outdir,outname1);
		//for (int i=0;i<input->size();i++)
		//{
		//	if (input->getXYZ(i)!=200)
		//	{
		//		input->putXYZ(i,0);
		//	}
		//}
		test.writeImageName(*initial,outdir);
	}
	else if (method == 3)//input data: inner data,outer data,no circle data,region growing stop when
	{


		


		
	}
	delete[] indata;

}
Seed *findplanepoint(Seed *outerPoint,Seed * centerPoint,int height)
{
	
	float gradient = sqrt((outerPoint->x - centerPoint->x) + (outerPoint->y - centerPoint->y) + (outerPoint->z - centerPoint->z));
	float deltax = float(outerPoint->x - centerPoint->x) / gradient;
	float deltay = float(outerPoint->y - centerPoint->y) / gradient;
	float deltaz = float(outerPoint->z - centerPoint->z) / gradient;
	int x = centerPoint->x + deltax*height;
	int y = centerPoint->y+deltay*height;
	int z = centerPoint->z + deltaz*height;
	Seed *planePoint = new Seed(x,y,z);
	return planePoint;
}
void findupplane(Raw *innerandouter,Seed *planePoint,Seed *normaldirection)
{
	vector<Seed> seedlist;
	int x = planePoint->x;
	int y = planePoint->y;
	int z = planePoint->z;
	for (size_t i = 0; i < innerandouter->getXsize(); i++)
	{
		for (size_t j = 0; j < innerandouter->getYsize(); j++)
		{
			for (size_t k = 0; k < innerandouter->getZsize(); k++)
			{
				if (innerandouter->get(i,j,k)!=0)
				{
					if ((i-x)*normaldirection->x+(j-y)*normaldirection->y+(k-z)*normaldirection->z==0)
					{
						innerandouter->put(i, j, k, 200);
					}
				}
			}
		}
	}
	//return innerandouter;
}
void geometrymethod()
{
	int l = 512, m = 512, n = 50;
	vector<Seed>seedlist;
	RawImage *test = new RawImage();
	int offset = 200;
	Raw *outer, *inner;
	PIXTYPE *outerdata = new PIXTYPE[l*m*n + offset];
	PIXTYPE *innerdata = new PIXTYPE[l*m*n + offset];
	outerdata = test->readStreamfloat("D:\\swfdata20140420res\\polyp\\3036P\\inner\\", &l, &m, &n);//outer wall 
	innerdata = test->readStreamfloat("D:\\swfdata20140420res\\polyp\\3036P\\outer\\", &l, &m, &n);//inner wall
	n = 50;
	outer = new Raw(l, m, n, outerdata);
	inner = new Raw(l, m, n, innerdata);
	Raw innerAndOuter = *outer * 2 + *inner;//boundary data,outer is 2,inner is 1
	ThreeDim_LevelSet *ls = new ThreeDim_LevelSet();
	ls->initialg(*outer);
	ls->initialg(*inner);
	Seed * middle = new Seed(169, 237, 11);
	seedlist.clear();
	seedlist.push_back(*middle);
	//next step draw a ball
	vector<Seed> outerlist = probabiltyouterpointleastcost(outer, middle, 10);
	Seed *leastcostseed = findLeastCostSeed(outerlist, middle);
	Seed *planePoint = findplanepoint(leastcostseed, middle, 8);//the height is the size of the polyp
	Seed *normal = new Seed(middle->x-planePoint->x,middle->y-planePoint->y,middle->z-planePoint->z);
	findupplane(&innerAndOuter,planePoint,normal);
	vector<Seed> start;
	Seed *seed = new Seed(169,237,11);
	start.push_back(*seed);
	seedGrowing(&innerAndOuter,start);




	
}
void seedGrowing(Raw * src,vector<Seed> seedlist)
{
	queue<Seed> seedqueue;
	//Seed seedl=;
	seedqueue.push(seedlist.back());
	seedlist.pop_back();

	while (!seedqueue.empty())
	{
		
		Seed start=seedqueue.front();
		seedqueue.pop();
		//int size=1;
		//for (int i=start.x-1; i <start.x+1;i++)
		//{
		//	for (int j=start.y-1; j<start.y+1;j++)
		//	{
		//		for (int k=start.z-1; k<start.z+1; k++)
		//		{
		//			if (src->get(i,j,k)==1)
		//			{
		//				Seed *adj=new Seed(i,j,k);
		//				seedqueue.push(*adj);
		//				src->put(i,j,k,200);
		//			}
		//		}
		//	}
		//}
		int i=start.x;
		int j=start.y;
		int k=start.z;
		//for (int i = 0; i < src->size(); i++)
		//{
		//	src->putXYZ(i,0);

		//}
		if (i>1&&i<src->getXsize()-1&&j>1 && j<src->getYsize()-1 && k > 1&& k<src->getZsize()-1)
		{
			int up=src->get(i-1,j,k)==1;
			int down=src->get(i+1,j,k)==1;
			int left=src->get(i,j-1,k)==1;
			int right=src->get(i,j+1,k)==1;
			int front =src->get(i,j,k-1)==1;
			int back=src->get(i,j,k+1)==1;
			if (src->get(i,j,k)==1)
			{
				int count=0;
				if (up)
				{
					Seed *adj=new Seed(i-1,j,k);
					seedqueue.push(*adj);
					//src->put(i-1,j,k,200);
					count++;
				}
				if (down)
				{
					Seed *adj=new Seed(i+1,j,k);
					seedqueue.push(*adj);
						count++;
					//src->put(i+1,j,k,200);
				}
				if (left)
				{
					Seed *adj=new Seed(i,j-1,k);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j-1,k,200);
				}
				if (right)
				{
					Seed *adj=new Seed(i,j+1,k);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j+1,k,200);
				}
				if (front)
				{
					Seed *adj=new Seed(i,j,k-1);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j,k-1,200);
				}
				if (back)
				{
					Seed *adj=new Seed(i,j,k+1);
					seedqueue.push(*adj);
						count++;
					//src->put(i,j,k+1,200);
				}
				//if (count>0)
				//{
					src->put(i,j,k,200);
				//}
				//else src->put(i,j,k,300);
			}
		}


	}//...while


}
void seedGrowing(Raw * src, vector<Seed> seedlist,int threshold)//threshold is the value to be put into the region
{
	queue<Seed> seedqueue;
	//Seed seedl=;
	seedqueue.push(seedlist.back());
	seedlist.pop_back();
	int value = threshold;
	while (!seedqueue.empty())
	{

		Seed start = seedqueue.front();
		seedqueue.pop();

		int i = start.x;
		int j = start.y;
		int k = start.z;

		if (i>1 && i<src->getXsize() - 1 && j>1 && j<src->getYsize() - 1 && k > 1 && k<src->getZsize() - 1)
		{
			int up = src->get(i - 1, j, k) == value;
			int down = src->get(i + 1, j, k) == value;
			int left = src->get(i, j - 1, k) == value;
			int right = src->get(i, j + 1, k) == value;
			int front = src->get(i, j, k - 1) == value;
			int back = src->get(i, j, k + 1) == value;
			if (src->get(i, j, k) == value)
			{
				int count = 0;
				if (up)
				{
					Seed *adj = new Seed(i - 1, j, k);
					seedqueue.push(*adj);
					//src->put(i-1,j,k,200);
					count++;
				}
				if (down)
				{
					Seed *adj = new Seed(i + 1, j, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i+1,j,k,200);
				}
				if (left)
				{
					Seed *adj = new Seed(i, j - 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j-1,k,200);
				}
				if (right)
				{
					Seed *adj = new Seed(i, j + 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j+1,k,200);
				}
				if (front)
				{
					Seed *adj = new Seed(i, j, k - 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k-1,200);
				}
				if (back)
				{
					Seed *adj = new Seed(i, j, k + 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k+1,200);
				}
				//if (count>0)
				//{
				src->put(i, j , k, 200);
				//}
				//else src->put(i,j,k,300);
			}
		}


	}//...while


}
void seedGrowing(Raw * src,Raw *origin, vector<Seed> seedlist)
{
	queue<Seed> seedqueue;
	//Seed seedl=;
	seedqueue.push(seedlist.back());
	seedlist.pop_back();

	while (!seedqueue.empty())
	{

		Seed start = seedqueue.front();
		seedqueue.pop();
		//int size=1;
		//for (int i=start.x-1; i <start.x+1;i++)
		//{
		//	for (int j=start.y-1; j<start.y+1;j++)
		//	{
		//		for (int k=start.z-1; k<start.z+1; k++)
		//		{
		//			if (src->get(i,j,k)==1)
		//			{
		//				Seed *adj=new Seed(i,j,k);
		//				seedqueue.push(*adj);
		//				src->put(i,j,k,200);
		//			}
		//		}
		//	}
		//}
		int i = start.x;
		int j = start.y;
		int k = start.z;
		int k2 = start.z + 70;
		//for (int i = 0; i < src->size(); i++)
		//{
		//	src->putXYZ(i,0);

		//}
		if (i>1 && i<src->getXsize() - 1 && j>1 && j<src->getYsize() - 1 && k > 1 && k<src->getZsize() - 1)
		{
			int value = 200;
			int originflag = origin->get(i,j,k)>-100;
			int up = src->get(i - 1, j, k) != value && origin->get(i-1, j, k2)>-100;
			//int val = origin->get(i - 1, j, k);
			int down = src->get(i + 1, j, k) != value && origin->get(i + 1, j, k2)>-100;
			int left = src->get(i, j - 1, k) != value && origin->get(i, j - 1, k2)>-100;
			int right = src->get(i, j + 1, k) != value && origin->get(i, j + 1, k2)>-100;
			int front = src->get(i, j, k - 1) != value && origin->get(i, j, k2 - 1)>-100;
			int back = src->get(i, j, k + 1) != value && origin->get(i, j, k2 + 1)>-100;
			if (src->get(i, j, k) != 200)
			{
				int count = 0;
				if (up)
				{
					Seed *adj = new Seed(i - 1, j, k);
					seedqueue.push(*adj);
					//src->put(i-1,j,k,200);
					count++;
				}
				if (down)
				{
					Seed *adj = new Seed(i + 1, j, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i+1,j,k,200);
				}
				if (left)
				{
					Seed *adj = new Seed(i, j - 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j-1,k,200);
				}
				if (right)
				{
					Seed *adj = new Seed(i, j + 1, k);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j+1,k,200);
				}
				if (front)
				{
					Seed *adj = new Seed(i, j, k - 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k-1,200);
				}
				if (back)
				{
					Seed *adj = new Seed(i, j, k + 1);
					seedqueue.push(*adj);
					count++;
					//src->put(i,j,k+1,200);
				}
				//if (count>0)
				//{
					src->put(i, j, k, 200);
				//}
				//else src->put(i, j, k, 300);
			}
		}


	}//...while


}