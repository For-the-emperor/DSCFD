
#include "stdafx.h"
#include <iostream>
#include "math.h"
#include <fstream>
#include <stdio.h>


using namespace std;


//�ó�������������
FILE *fp=NULL;    //STL�ļ������룩
char fpName[200]="D:\\Documents\\Electromagnetic\\model\\��\\Simplified Airplane\\airplane.stl";  //STL�ļ�λ��
FILE *fp2=NULL;   //txt�ļ���д����
char fp2Name[200]="D:\\Documents\\Electromagnetic\\model\\��\\Simplified Airplane\\airplane_1.5G.txt";   //����д����ļ�λ��
int counttol=0;   //Ŀ����Ԫ���������-1
double *xres;      //Ŀ����Ԫ����������
double grid_minx; double grid_maxx;
double grid_miny; double grid_maxy;
double grid_minz; double grid_maxz;

//Ԫ����Ϣ
double c_0 = 2.997924580105029e+08;
int N_per_wavelength = 15;
double f = 1.5e+09;
double lemda = c_0 / f;
double dx = lemda / N_per_wavelength;
double dy = lemda / N_per_wavelength;
double dz = lemda / N_per_wavelength;//Ԫ���ߴ�

//������Ƭ��Ϣ
int n_trigon=0;//����
double *x;
double *y;
double *z;//��������
double *cout1;//����������
double *pcor;

//��������
double min(double *arry,int size);
double max(double *arry,int size);//һά������ֵ
double min2s(double a,double b);
double max2s(double a,double b);//����������ֵ
int xiabiao(double *x,int size,double y);//abs(x[]-y)����Сֵ���±�
int any(double *a);//�ж�һά3Ԫ�������Ƿ��з�0Ԫ�أ��з���1��û�з���0
void cross(double *A,double*B);//������ά���������飩�������������
double sum(double *A,double*B);//����3Ԫ�������ӦԪ���������ͣ�������(��)����
double normest(double *A,int size);//һά������2����
double *chu(double *p,int size,double a);//����p���Գ���a
double *sort(double *p,int size);//������������

void read();//��stlASCLL�ļ�
void poufen();//�ʷ�
void write();//��xresд��txt

//������
int main()
{
	read();
	printf("������Ƭ������%d\n",n_trigon);
	poufen();
	printf("Ŀ����Ԫ�����������counttol=%d\n",counttol);
	write();
	printf("����ļ�λ�ã�%s\n",fp2Name);
	return 0;
}

//��STL�ļ�����
void read()
{
	char partName[80];//�ļ�ͷ
 
	struct Normal  //����������
	{
		double i;
		double j;
		double k;
	};
	struct Vertrex  //��������
	{
		double x;
		double y;
		double z;
	};
	struct Face
	{
		Normal normal; //������
		Vertrex vertex1; //����
		Vertrex vertex2;
		Vertrex vertex3;
	};
 
	//��������Ƭ���������
	struct Faces
	{ 
		Face data;
		Faces *pointer;  //ָ����ָ��ڵ��ָ�룩
	};
	
	//FILE *fp=NULL;
	Faces *faces = (Faces*)malloc(sizeof(Faces));//����ı�ͷ
	Faces *p = faces; //p�ȼ���faces,ָ��ṹ��Faces�ı������׵�ַ
	Faces *q = NULL;
	fp = fopen(fpName, "r"); //���ļ�·�����ļ�
	//printf("fp=%d\n",fp);
 
	if (fp)
	{
		char str[60];
		fscanf(fp, "%s", str);//�ô�������ý����ƶ��ļ�ָ�룬ʹָ��������ʶ����solid�������ļ�ָ���Ƶ��ļ����Ƶ�λ��
		fscanf(fp, "%s", partName);//�����ļ�����
		printf("�ļ�ͷ��%s\n",partName);

		//ѭ����ȡ����Ƭ��ķ��ߺͶ�������
		while (!feof(fp)) //�ļ����ݶ���֮ǰ��feof��ֵ��Ϊ0
		{
			p->pointer= (Faces*)malloc(sizeof(Faces));//�����½�㣬���ڽ�������  pointer��p����һ����ַ
			q = p;
			p = p->pointer;
			fscanf(fp, "%s%s", str,str);//ʹ�ļ�ָ��������ʶ����facet normal����ָ���������ݲ���
			fscanf(fp, "%lf%lf%lf", &p->data.normal.i, &p->data.normal.j, &p->data.normal.k);//��ȡ��������
			fscanf(fp, "%s%s", str,str);//ʹ�ļ�ָ��������ʶ����outer loop��
			fscanf(fp, "%s", str);//ʹ�ļ�ָ��������ʶ����vertex����ָ������1�����ݲ���
			fscanf(fp, "%lf%lf%lf", &p->data.vertex1.x, &p->data.vertex1.y, &p->data.vertex1.z);
			fscanf(fp, "%s", str);//ʹ�ļ�ָ��������ʶ����vertex����ָ������2�����ݲ���
			fscanf(fp, "%lf%lf%lf", &p->data.vertex2.x, &p->data.vertex2.y, &p->data.vertex2.z);
			fscanf(fp, "%s", str);//ʹ�ļ�ָ��������ʶ����vertex����ָ������3�����ݲ���
			fscanf(fp, "%lf%lf%lf", &p->data.vertex3.x, &p->data.vertex3.y, &p->data.vertex3.z);
			fscanf(fp, "%s", str);//ʹ�ļ�ָ��������ʶ����endloop��
			fscanf(fp, "%s", str);//ʹ�ļ�ָ��������ʶ����endfacet��
			n_trigon++;
			
		}
		//printf("%d\n%d\n",p,q->pointer); //p��q->pointer��ַһ����Ϊʲô���������д�����p����
		free(q->pointer);//�����ļ�ĩβ��һ���ַ���endsolid .....������ɶഴ����һ����㣬�ý������������õģ����Ҫ�ͷ�
		q->pointer = NULL;//���һ�����������ͷ���֮�󣬽�ָ��������ָ����ΪNULL�������������
		n_trigon=n_trigon-1;
	    //printf("������Ƭ������%d\n",n_trigon);
	}
	else
	{
		printf("���ļ�ʧ�ܣ�");
	}

    //Ϊ��������ͷ��������䶯̬�ڴ�
	cout1=(double*)malloc(sizeof(double)*n_trigon*3);
	x=(double*)malloc(sizeof(double)*n_trigon*3);
	y=(double*)malloc(sizeof(double)*n_trigon*3);
	z=(double*)malloc(sizeof(double)*n_trigon*3);
	p = faces->pointer;//��ȡ����ڶ�������ָ��
	int li=0;
	while (p!=NULL)
	{   //��ֵ����ʹ��
		cout1[li*3+0]=p->data.normal.i;
		cout1[li*3+1]=p->data.normal.j;
		cout1[li*3+2]=p->data.normal.k;

		x[li*3+0]=p->data.vertex1.x;
		x[li*3+1]=p->data.vertex2.x;
		x[li*3+2]=p->data.vertex3.x;

		y[li*3+0]=p->data.vertex1.y;
		y[li*3+1]=p->data.vertex2.y;
		y[li*3+2]=p->data.vertex3.y;

		z[li*3+0]=p->data.vertex1.z;
		z[li*3+1]=p->data.vertex2.z;
		z[li*3+2]=p->data.vertex3.z;
			
		//ѭ�����
		//printf("���ߣ�\n");
		//printf("%f  %f  %f\n", cout1[li*3+0], cout1[li*3+1], cout1[li*3+2]);
		//printf("���㣺\n");
		//printf("%f  %f  %f\n", x[li*3+0], y[li*3+0], z[li*3+0]);
		//printf("%f  %f  %f\n", x[li*3+1], y[li*3+1], z[li*3+1]);
		//printf("%f  %f  %f\n", x[li*3+2], y[li*3+2], z[li*3+2]);

		li++;
		p = p->pointer;
		//printf("%d\n",li);
	}
}

//�ʷֳ���
void poufen()
{
	double xmin=0;
	double xmax=0;
	double ymin=0;
	double ymax=0;
	double zmin=0;
	double zmax=0;

	int x_size=n_trigon*3;
	int y_size=n_trigon*3;
	int z_size=n_trigon*3;
	xmin=min(x,x_size);xmax=max(x,x_size);
	ymin=min(y,y_size);ymax=max(y,y_size);
	zmin=min(z,z_size);zmax=max(z,z_size);
	//printf("xmin:%f  ymin:%f  zmin:%f\n",xmin,ymin,zmin);
	//printf("xmax:%f  ymax:%f  zmax:%f\n",xmax,ymax,zmax);
	//for(int i=0;i<n_trigon*3;i++)
	//	{
	//		x[i]=x[i]-xmin;
	//		y[i]=y[i]-ymin;
	//		z[i]=z[i]-zmin;
	//	}
	//xmax=xmax-xmin;xmin=0;
	//ymax=ymax-ymin;ymin=0;
	//zmax=zmax-zmin;zmin=0;
	//xmin=min(x,x_size);xmax=max(x,x_size);
	//ymin=min(y,y_size);ymax=max(y,y_size);
	//zmin=min(z,z_size);zmax=max(z,z_size);
	grid_minx = ((int)(xmin / dx))*dx; grid_maxx = ((int)(xmax / dx))*dx;
	grid_miny = ((int)(ymin / dy))*dy; grid_maxy = ((int)(ymax / dy))*dy;
	grid_minz = ((int)(zmin / dz))*dz; grid_maxz = ((int)(zmax / dz))*dz;
	printf("grid_minx:%f  grid_miny:%f  grid_minz:%f\n", grid_minx, grid_miny, grid_minz);
	printf("grid_maxx:%f  grid_maxy:%f  grid_maxz:%f\n", grid_maxx, grid_maxy, grid_maxz);

	int n_xx = (grid_maxx - grid_minx) / dx + 1;
	int n_yy = (grid_maxy - grid_miny) / dy + 1;
	int n_zz = (grid_maxz - grid_minz) / dz + 1;
	printf("3�������Ԫ��������nxx=%d  nyy=%d  nzz=%d\n",n_xx,n_yy,n_zz);

	double *xx=(double *)malloc(sizeof(double)*(n_xx));
	double *yy=(double *)malloc(sizeof(double)*(n_yy));
	double *zz=(double *)malloc(sizeof(double)*(n_zz));//Ԫ����������
	for (int i = 0; i<n_xx; i++)
	{
		xx[i] = grid_minx + i*dx;
	}
	for (int i = 0; i<n_yy; i++)
	{
		yy[i] = grid_miny + i*dy;
	}
	for (int i = 0; i<n_zz; i++)
	{
		zz[i] = grid_minz + i*dz;
	}
	

	//��ʼ��
	double x2s[2]={0,0};
	double xmin_2s=0;
	double xmax_2s=0;
	double xf12min=0;
	double xf12max=0;

	double jdu=0;//������Ƭ�����������ߴ�ֱʱΪ0
	int jd=0;//�������ޣ���Ϊ1����Ϊ0
	int as1=0;
	int as2=0;//�����±�
	int count=0;//�������
	int count1=0;
	counttol=0;//Ŀ����Ԫ���������
	int zeroflag=0;
	double xf[3]={0,0,0};
	double diff=1e-10;

	double A1=0;double B1=0;double C1=0;double xe=0;
	double AB[3]={0,0,0};double AC[3]={0,0,0};double AP[3]={0,0,0};
	double BA[3]={0,0,0};double BC[3]={0,0,0};double BP[3]={0,0,0};
	double CA[3]={0,0,0};double CB[3]={0,0,0};double CP[3]={0,0,0};

	double *jk1;double *jk2;double *jk3;double *jk4;double *jk5;double *jk6;
	jk1=(double *)malloc(sizeof(double)*3);jk2=(double *)malloc(sizeof(double)*3);
	jk3=(double *)malloc(sizeof(double)*3);jk4=(double *)malloc(sizeof(double)*3);
	jk5=(double *)malloc(sizeof(double)*3);jk6=(double *)malloc(sizeof(double)*3);
	pcor=(double *)malloc(sizeof(double)*3);
	for(int i=0;i<3;i++)
	{
		jk1[i]=0;jk2[i]=0;jk3[i]=0;jk4[i]=0;jk5[i]=0;jk6[i]=0;
	}

	double *xd; //λ��Ŀ�������ڲ�����������Ƭ���������Ԫ��x����
	xd=(double*)malloc(sizeof(double)*(n_xx));
    for(int jm=0;jm<n_xx;jm++)
	{
		xd[jm]=0;
	}
	int len_xd=0;//xd��ЧԪ�ظ���
	xres=(double*)malloc(sizeof(double)*n_xx*n_yy*n_zz*3);//���Ŀ����Ԫ����������
	for(int i=0;i<n_xx*n_yy*n_zz*3;i++)
	{
		xres[i]=0;
	}

	double x1=0;
	//��ʼѭ��
	//����x��ƽ�е�����������������Ƭ�Ľ���
	int n_ver=n_trigon;
	for(int i=0;i<n_yy;i++)
	{
		for(int j=0;j<n_zz;j++)
		{
			jd=0;
			//count=0;//���µ����ߣ�count��Ҫ��ʼ��Ϊ0
			for(int m=0;m<n_ver;m++)
			{
				jdu=cout1[m*3+0]*(-x[m*3+0])+cout1[m*3+1]*(yy[i]-y[m*3+0])+cout1[m*3+2]*(zz[j]-z[m*3+0]);
				if((fabs(cout1[m*3+0])<diff) && (fabs(jdu)<diff))//��ʱ����������
				{
					//��12���غ�
					if(fabs(y[m*3+0]-y[m*3+1])<diff && (fabs(y[m*3+0]-yy[i]))<diff && (fabs(z[m*3+0]-z[m*3+1]))<diff && (fabs(z[m*3+0]-zz[j]))<diff)
					{
						x2s[0]=x[m*3+0];
						x2s[1]=x[m*3+1];
						xmin_2s=min(x2s,2);
						xmax_2s=max(x2s,2);
						as1=xiabiao(xx,n_xx,xmin_2s);
						as2=xiabiao(xx,n_xx,xmax_2s);
						jd=1;
						for(int k=0;k<(as2-as1+1);k++)
						{
							xd[count]=xx[as1]+dx*k;
							count++;
						}
						zeroflag=1;
					}
					//��23���غ�
					if(fabs(y[m*3+1]-y[m*3+2])<diff && (fabs(y[m*3+1]-yy[i]))<diff && (fabs(z[m*3+1]-z[m*3+2]))<diff && (fabs(z[m*3+1]-zz[j]))<diff)
					{
						x2s[0]=x[m*3+1];
						x2s[1]=x[m*3+2];
						xmin_2s=min(x2s,2);
						xmax_2s=max(x2s,2);
						as1=xiabiao(xx,n_xx,xmin_2s);
						as2=xiabiao(xx,n_xx,xmax_2s);
						jd=1;
						for(int k=0;k<(as2-as1+1);k++)
						{
							xd[count]=xx[as1]+dx*k;
							count++;
						}
						zeroflag=1;
					}
					//��13���غ�
					if(fabs(y[m*3+0]-y[m*3+2])<diff && (fabs(y[m*3+0]-yy[i]))<diff && (fabs(z[m*3+0]-z[m*3+2]))<diff && (fabs(z[m*3+0]-zz[j]))<diff)
					{
						x2s[0]=x[m*3+0];
						x2s[1]=x[m*3+2];
						xmin_2s=min(x2s,2);
						xmax_2s=max(x2s,2);
						as1=xiabiao(xx,n_xx,xmin_2s);
						as2=xiabiao(xx,n_xx,xmax_2s);
						jd=1;
						for(int k=0;k<(as2-as1+1);k++)
						{
							xd[count]=xx[as1]+dx*k;
							count++;
						}
						zeroflag=1;
					}
					if(zeroflag==0)//�����߶����غ�
					{
						//float x1=0;
						//1��
						if((fabs(y[m*3+1]-y[m*3+0])<diff)&&(fabs(z[m*3+1]-z[m*3+0])>diff))
						{
							x1=(zz[j]-z[m*3+0])*(x[m*3+1]-x[m*3+0])/(z[m*3+1]-z[m*3+0])+x[m*3+0];
							if(x1>=min2s(x[m*3+0],x[m*3+1])&&x1<=max2s(x[m*3+0],x[m*3+1]))
							{
								jd=1;
								xf[count1]=x1;
								count1++;
								xf[2]=1;
							}
						}
						if(fabs(y[m*3+1]-y[m*3+0])>diff)
						{
							x1=(yy[i]-y[m*3+0])*(x[m*3+1]-x[m*3+0])/(y[m*3+1]-y[m*3+0])+x[m*3+0];
							if(x1>=min2s(x[m*3+0],x[m*3+1])&&x1<=max2s(x[m*3+0],x[m*3+1]))
							{
								jd=1;
								xf[count1]=x1;
								count1++;
								xf[2]=1;
							}
						}
						//2��
						if((fabs(y[m*3+2]-y[m*3+0])<diff)&&(fabs(z[m*3+2]-z[m*3+0])>diff))
						{
							x1=(zz[j]-z[m*3+0])*(x[m*3+2]-x[m*3+0])/(z[m*3+2]-z[m*3+0])+x[m*3+0];
							if(x1>=min2s(x[m*3+0],x[m*3+2])&&x1<=max2s(x[m*3+0],x[m*3+2]))
							{
								jd=1;
								xf[count1]=x1;
								count1++;
								xf[2]=1;
							}
						}
						if(fabs(y[m*3+2]-y[m*3+0])>diff)//2��
						{
							x1=(yy[i]-y[m*3+0])*(x[m*3+2]-x[m*3+0])/(y[m*3+2]-y[m*3+0])+x[m*3+0];
							if(x1>=min2s(x[m*3+0],x[m*3+2])&&x1<=max2s(x[m*3+0],x[m*3+2]))
							{
								jd=1;
								xf[count1]=x1;
								count1++;
								xf[2]=1;
							}
						}
						//3��
						if((fabs(y[m*3+2]-y[m*3+1])<diff)&&(fabs(z[m*3+2]-z[m*3+1])>diff))
						{
							x1=(zz[j]-z[m*3+1])*(x[m*3+2]-x[m*3+1])/(z[m*3+2]-z[m*3+1])+x[m*3+1];
							if(x1>=min2s(x[m*3+1],x[m*3+2])&&x1<=max2s(x[m*3+1],x[m*3+2]))
							{
								jd=1;
								xf[count1]=x1;
								count1++;
								xf[2]=1;
							}
						}
						if(fabs(y[m*3+2]-y[m*3+1])>diff)//3��
						{
							x1=(yy[i]-y[m*3+1])*(x[m*3+2]-x[m*3+1])/(y[m*3+2]-y[m*3+1])+x[m*3+1];
							if(x1>=min2s(x[m*3+1],x[m*3+2])&&x1<=max2s(x[m*3+1],x[m*3+2]))
							{
								jd=1;
								xf[count1]=x1;
								count1++;
								xf[2]=1;
							}
						}

						if(xf[2]==1 && count1==2)
						{
							jd=1;
							xf12min=min2s(xf[0],xf[1]);
							xf12max=max2s(xf[0],xf[1]);
							as1=xiabiao(xx,n_xx,xf12min);
							as2=xiabiao(xx,n_xx,xf12max);
							for(int k=0;k<(as2-as1+1);k++)
							{
								xd[count]=xx[as1]+dx*k;
								count++;
							}
						}
						count1=0;
						for(int k=0;k<3;k++)
						{
							xf[k]=0;
						}
					}
					zeroflag=0;
				}
				else//��ʱ�߲������ϣ�ֻ��һ�������û�н���
				{
					A1=(y[m*3+1]-y[m*3+0])*(z[m*3+2]-z[m*3+0])-(y[m*3+2]-y[m*3+0])*(z[m*3+1]-z[m*3+0]);
					if(fabs(A1)>1e-8)
					{
						B1=-(x[m*3+1]-x[m*3+0])*(z[m*3+2]-z[m*3+0])+(x[m*3+2]-x[m*3+0])*(z[m*3+1]-z[m*3+0]);
						C1=(x[m*3+1]-x[m*3+0])*(y[m*3+2]-y[m*3+0])-(x[m*3+2]-x[m*3+0])*(y[m*3+1]-y[m*3+0]);
						xe=(-B1*(yy[i]-y[m*3+0])-C1*(zz[j]-z[m*3+0]))/A1+x[m*3+0];

						AC[0]=x[m*3+2]-x[m*3+0];AC[1]=y[m*3+2]-y[m*3+0];AC[2]=z[m*3+2]-z[m*3+0];
						AP[0]=xe-x[m*3+0];		AP[1]=yy[i]-y[m*3+0];	AP[2]=zz[j]-z[m*3+0];
						AB[0]=x[m*3+1]-x[m*3+0];AB[1]=y[m*3+1]-y[m*3+0];AB[2]=z[m*3+1]-z[m*3+0];

						BA[0]=x[m*3+0]-x[m*3+1];BA[1]=y[m*3+0]-y[m*3+1];BA[2]=z[m*3+0]-z[m*3+1];
						BP[0]=xe-x[m*3+1];		BP[1]=yy[i]-y[m*3+1];	BP[2]=zz[j]-z[m*3+1];
						BC[0]=x[m*3+2]-x[m*3+1];BC[1]=y[m*3+2]-y[m*3+1];BC[2]=z[m*3+2]-z[m*3+1];

						CB[0]=x[m*3+1]-x[m*3+2];CB[1]=y[m*3+1]-y[m*3+2];CB[2]=z[m*3+1]-z[m*3+2];
						CP[0]=xe-x[m*3+2];		CP[1]=yy[i]-y[m*3+2];	CP[2]=zz[j]-z[m*3+2];
						CA[0]=x[m*3+0]-x[m*3+2];CA[1]=y[m*3+0]-y[m*3+2];CA[2]=z[m*3+0]-z[m*3+2];

						cross(AP,AB);
						jk1[0] = pcor[0]; jk1[1] = pcor[1]; jk1[2] = pcor[2]; 
						cross(AC,AB);
						jk2[0] = pcor[0]; jk2[1] = pcor[1]; jk2[2] = pcor[2]; 
						cross(BP,BC);
						jk3[0] = pcor[0]; jk3[1] = pcor[1]; jk3[2] = pcor[2]; 
						cross(BA,BC);
						jk4[0] = pcor[0]; jk4[1] = pcor[1]; jk4[2] = pcor[2]; 
						cross(CP,CA);
						jk5[0] = pcor[0]; jk5[1] = pcor[1]; jk5[2] = pcor[2]; 
						cross(CB,CA);
						jk6[0] = pcor[0]; jk6[1] = pcor[1]; jk6[2] = pcor[2]; 
						
						if((((!any(jk1))&&(sum(AP,AB)>=0))&&(normest(AP,3)<normest(AB,3)))||
							(((!any(jk3))&&(sum(BP,BC)>=0))&&(normest(BP,3)<normest(BC,3)))||
							(((!any(jk5))&&(sum(CP,CA)>=0))&&(normest(CP,3)<normest(CA,3))))
						{
							jd=1;
							as1=xiabiao(xx,n_xx,xe);
							xd[count]=xx[as1];
							count++;
						}
						else
						{
							jk1=chu(jk1,3,sqrt(sum(jk1,jk1)));jk2=chu(jk2,3,sqrt(sum(jk2,jk2)));
							jk3=chu(jk3,3,sqrt(sum(jk3,jk3)));jk4=chu(jk4,3,sqrt(sum(jk4,jk4)));
							jk5=chu(jk5,3,sqrt(sum(jk5,jk5)));jk6=chu(jk6,3,sqrt(sum(jk6,jk6)));
							if((sum(jk1,jk2)>=0.5)&&(sum(jk3,jk4)>=0.5)&&(sum(jk5,jk6)>=0.5))
							{
								jd=1;
								as1=xiabiao(xx,n_xx,xe);
								xd[count]=xx[as1];
								count++;
							}
						}	
					}
				}	
			}
			
			if(jd==1)
			{
				len_xd=count;
				xd=sort(xd,len_xd);//��������
				if(len_xd==1)
				{
					xres[counttol*3+0]=xd[0];
					xres[counttol*3+1]=yy[i];
					xres[counttol*3+2]=zz[j];
					counttol++;
				}
				else
				{
					if(len_xd%2==0)//ż��
					for(int kk=0;kk<len_xd;)
					{
						int NN=(int)((xd[kk+1]-xd[kk])/dx+1+0.5);//��������
						for(int MM=counttol;MM<counttol+NN;MM++)
						{
							xres[MM*3+0]=xd[kk]+dx*(MM-counttol);
							xres[MM*3+1]=yy[i];
							xres[MM*3+2]=zz[j]; 
						}
						counttol=counttol+NN;
						kk=kk+2;
					}
					else//����
					{
						for(int kk=0;kk<len_xd-1;)
						{
							int NN=(int)((xd[kk+1]-xd[kk])/dx+1+0.5);
							for(int MM=counttol;MM<counttol+NN;MM++)
							{
								xres[MM*3+0]=xd[kk]+dx*(MM-counttol);
								xres[MM*3+1]=yy[i];
								xres[MM*3+2]=zz[j];
							}
							counttol=counttol+NN;
							kk=kk+2;
						}
						xres[counttol*3+0]=xd[len_xd-1];
						xres[counttol*3+1]=yy[i];
						xres[counttol*3+2]=zz[j];
						counttol++;
					}
				}
				for(int i=0;i<len_xd;i++)
				{
					xd[i]=0;
				}
				count=0;
			}
			printf("i=%d of %d,j=%d of %d  counttol=%d",i,n_yy,j,n_zz,counttol);
			printf("  len_xd=%d\n",len_xd);
		}	
	}
	//printf("counttol=%d\n",counttol);
	//�������
	printf("Ŀ�꣺\n");
	printf("jdu=%f\n",jdu);
	printf("count=%d  count1=%d  counttol=%d  zeroflag=%d\n",count,count1,counttol,zeroflag);

    /*printf("A1=%f   B1=%f   C1=%f   xe=%f\n",A1,B1,C1,xe);
	for(int i=0;i<3;i++)
	{
		printf("AB[%d]=%f  AC[%d]=%f  AP[%d]=%f\n",i,AB[i],i,AC[i],i,AP[i]);
	}
	for(int i=0;i<3;i++)
	{
		printf("BA[i]=%f  BC[i]=%f  BP[i]=%f\n",BA[i],BC[i],BP[i]);

	}
	for(int i=0;i<3;i++)
	{
		printf("CA[i]=%f  CB[i]=%f  CP[i]=%f\n",CA[i],CB[i],CP[i]);
	}

	for(int i=0;i<3;i++)
	{
		printf("jk1[i]=%f  jk2[i]=%f  jk3[i]=%f\n",jk1[i],jk2[i],jk3[i]);
	}
	for(int i=0;i<3;i++)
	{
		printf("jk4[i]=%f  jk5[i]=%f  jk6[i]=%f\n",jk4[i],jk5[i],jk6[i]);
	}

	printf("len_xd=%d\n",len_xd);
	for(int i=0;i<len_xd;i++)
	{
		printf("xd[%d]=%f\n",i,xd[i]);
	}*/
}

//д��txt�ļ�����
void write()
	{
		//FILE *fp2=NULL;
		fp2=fopen(fp2Name, "w");
		if(fp2==NULL)
		{
			printf("Can not open txt file!\n");
			exit(1);
		}
		for(int i=0;i<counttol;i++)
		{
			fprintf(fp2,"%.10f   %.10f   %.10f\n",xres[i*3+0], xres[i*3+1], xres[i*3+2]);
		}
		fclose(fp2);
}

//С����

//��һά������ֵ
double min(double *arry,int size)
{
	double min=arry[0];
	int i;
	for(i=0;i<size;i++)
			if(min>arry[i])
				min=arry[i];
	return min;
}
double max(double *arry,int size)
{
	double max=arry[0];
	int i;
	for(i=0;i<size;i++)
			if(max<arry[i])
				max=arry[i];
	return max;
}

//��������ֵ
double min2s(double a,double b)
{
	double min=a;
	if(a>b)
		min=b;
	return min;
}
double max2s(double a,double b)
{
	double max=a;
	if(a<b)
		max=b;
	return max;
}

//��min(fabs(x-y))���±꣬x�����飬size��x�Ĵ�С��y�ǳ���
int xiabiao(double *x,int size,double y)
{
	double *x1=(double *)malloc(sizeof(double)*size);
	for(int i=0;i<size;i++)
		{
			x1[i]=fabs(x[i]-y);
		}
	int as=0;
	double min_abs=x1[0];
	for(int i=0;i<size;i++)
	{
		if(x1[i]<min_abs)
		{
			min_abs=x1[i];
			as=i;
		}	
	}
	return as;
}

//any�������ж�һά3Ԫ�������Ƿ��з�0Ԫ�أ��з���1��û�з���0��
int any(double *a)
 {
	 int result=0;
	 for(int i=0;i<3;i++)
	 {
		 if(a[i]!=0)
		 {
			 result =1;
			 break;
		 }
			 
	 }
	 return result;
 }

//����3Ԫ�������ӦԪ���������ͣ�������(��)����
double sum(double *A,double*B)
{
	double SUM=0;
	for(int i=0;i<3;i++)
		SUM=SUM+A[i]*B[i];
	return SUM;
}

//��������λ�����Ĳ��
void cross(double *A,double *B)
{
	pcor[0]=A[1]*B[2]-A[2]*B[1];
	pcor[1]=-A[0]*B[2]+A[2]*B[0];
	pcor[2]=A[0]*B[1]-A[1]*B[0];
}

//��һά������2����
double normest(double *A,int size)
{
	double result=0;
	for(int i=0;i<size;i++)
	{
		result+=fabs(A[i])*fabs(A[i]);
	}
	result=sqrt(result);
	return result;
}

//������Գ��� p[size]/a
double *chu(double *p,int size,double a)
{
	for(int i=0;i<size;i++)
	{
		p[i]=p[i]/a;
	}
	return p;
}

//����ǰsize��Ԫ����������
double *sort(double *p,int size)
{
	for(int j=0;j<size-1;j++)
	{
		for(int i=j+1;i<size;i++)
		{
			if(p[i]<p[j])
			{
				double a;
				a=p[i];
				p[i]=p[j];
				p[j]=a;
			}
		}
	}
	return p;
}
