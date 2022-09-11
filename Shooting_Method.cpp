#include<bits/stdc++.h>
#define g(x,y,y_pr) ((32+2*pow(x,3)-y*y_pr)/8) 
#define f(x,y,y_pr) (y_pr)
#define tol pow(10,-9)
#define ana(x) (x*x+(16/x))
#define itmax 50
using namespace std;

double p(double x,double psi,double psi_pr,double y,double y_pr)
{
	return ((-psi*y_pr-psi_pr*y)/8);
}
double q(double x,double psi,double psi_pr)
{
	return psi_pr;
}

double *rk4(double x,double y,double y_pr,double psi,double psi_pr,double h,double a[])
{
	double *y_psi;
	y_psi=new double[2];
	for(int i=0;i<20;i++)
	{
	double k1=0,k2=0,k3=0,k4=0,l1=0,l2=0,l3=0,l4=0,y0=0,z0=0;

	k1=h*f(x,y,y_pr);
	l1=h*g(x,y,y_pr);
	
	k2=h*f((x+h/2),(y+k1/2),(y_pr+l1/2));
	l2=h*g((x+h/2),(y+k1/2),(y_pr+l1/2));
	
	k3=h*f((x+h/2),(y+k2/2),(y_pr+l2/2));
	l3=h*g((x+h/2),(y+k2/2),(y_pr+l2/2));
	
	k4=h*f((x+h),(y+k3),(y_pr+l3));
	l4=h*g((x+h),(y+k3),(y_pr+l3));
	
	y0=(k1+2*k2+2*k3+k4)/6;
	z0=(l1+2*l2+2*l3+l4)/6;
	y+=y0;
    y_pr+=z0;
    

    
	k1=h*q(x,psi,psi_pr);
	l1=h*p(x,psi,psi_pr,y,y_pr);
	
	k2=h*q((x+h/2),(psi+k1/2),(psi_pr+l1/2));
	l2=h*p((x+h/2),(psi+k1/2),(psi_pr+l1/2),y,y_pr);
	
	k3=h*q((x+h/2),(psi+k2/2),(psi_pr+l2/2));
	l3=h*p((x+h/2),(psi+k2/2),(psi_pr+l2/2),y,y_pr);
	
	k4=h*q((x+h),(psi+k3),(psi_pr+l3));
	l4=h*p((x+h),(psi+k3),(psi_pr+l3),y,y_pr);
	
	y0=(k1+2*k2+2*k3+k4)/6;
	z0=(l1+2*l2+2*l3+l4)/6;
	psi+=y0;
    psi_pr+=z0;
    
    x+=h;
	a[i]=y;
	}
   y_psi[0]=y;
   y_psi[1]=psi;
   
   return y_psi;    

}

int main()
{
	double x,y,y_pr,psi,psi_pr,h=0.1,error=1,*y_psi,*yexact;
 	y_psi=new double[2];
 	yexact=new double[20];
	int count=0;
	ofstream output;
	output.open("assignment9.txt");
	output<<setprecision(5);


	cout<<"enter the initial guess of m : ";
	cin>>y_pr;//the value of m
	output<<endl<<endl<<"x\t"<<"Y\t\t"<<"Yana"<<endl;
		do
		{
			
		x=1.00,y=17.00,psi=0,psi_pr=1.00;
		y_psi=rk4(x,y,y_pr,psi,psi_pr,h,yexact);
			
		error=(y_psi[0]-(14.3333));	
	    y_pr -= (error/(y_psi[1]));//newton raphson
		count++;	
	
		}while(fabs(error)>tol );
	
	x=1.1;	
	for(int i=0;i<20;i++)
	{
		output<<x<<"\t"<<yexact[i]<<"\t"<<ana(x)<<endl;
		x+=0.1;
	}
		
		if(error>tol)
		{
			output<<"the solution didnot converge.";
		}
	
	output<<"The number of iterations it took:  "<<count<<endl;
	cout<<"Kindly check the output file";
	output.close();
	return 0;
}