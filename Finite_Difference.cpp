#include<bits/stdc++.h>
#define a(x,h) (1+(cos(x)*h/2))//assigning functions for element in 'a' matrix
#define d(x,h) (-2-sin(x)*h*h)
#define c(x,h) (1-(cos(x)*h/2))
#define b_vec(x,h) (-exp(x)*h*h)
#define n 10
using namespace std;
double y[n],a[n][n];

void TDMA(double x,double h,double ratio[])
{

int i,j,k;
double b[n],y0=0,yn=1;

for(i=1;i<=n-1;i++)//generating b vector
{
    if(i==1)
    {
    b[i]=b_vec(x,h)-a(x,h)*y0;
    }
    else if(i==n-1)
    {
    b[i]=b_vec(x,h)-c(x,h)*yn;
    }
    else
	b[i]=b_vec(x,h);

    x+=h;
}

for(k=2;k<=n-1;k++)//gauss elimination
{
	b[k]= b[k]- ratio[k]* b[k-1];
}

y[n-1]=b[n-1]/a[n-1][n-1];//back substitution
for(i=(n-2); i>=1; i--)
	y[i]=(b[i]-a[i][i+1]*y[i+1])/a[i][i];



}



int main()
{

int i,j,k;
double xn=1,x0=0,h=(xn-x0)/n,x=h,ratio[n];
h=(xn-x0)/n;
ofstream sol;
sol.open("solution_n_10.txt");

for(i=1;i<=n;i++)//generating A matrix
{
	if(i==1)//first row
	{
		a[i][i]=d(x,h);
		a[i][i+1]=c(x,h);
	}
	else if(i==n)//last row
	{
		a[i][i]=d(x,h);
		a[i][i-1]=a(x,h);
	}
	else//other rows
	{
		a[i][i]=d(x,h);
		a[i][i+1]=c(x,h);
		a[i][i-1]=a(x,h);
	}

	x+=h;
}

for(k=2;k<=n-1;k++)//gauss elimination
{
	ratio[k]=(a[k][k-1]/a[k-1][k-1]);
	a[k][k] = a[k][k] - ratio[k] * a[k-1][k];
	a[k][k-1]=0;
}
    x=h;

	TDMA(x,h,ratio);


	sol<<"y= [";
	for(i=1;i<=n-1;i++)//printing y vector
	{
		sol<<"\n"<<y[i]<<" ";
	}
	sol<<"]";
	sol<<endl<<endl;


 cout<<"Check the output file 'solution.txt'";
 return 0;
}
