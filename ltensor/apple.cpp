#include "apple.h"


apple::apple(int ww):
	vec(6,0.)
{
	weight = ww;
}
apple::apple():
	vec(6,0.)
{
	weight = 0. ;
}

apple::~apple()
{

}
	
int apple::setWeight(double wei){
	weight = wei ;
	double cc{1.};
	vec[0] = cc;


	vector<double> ret(6,0.);
	vector<vector<double> > matrix(6,vector<double>(6,0.)); 

	matrix[0][0] =123;
	matrix[2][4] =23;

	for (int i = 0; i < 6; ++i)
	{
		for (int j = 0; j < 6; ++j)
		{
			ret[i] += matrix[i][j]* vec[j];
		}
	}

	cout<<"ret[0] "<<ret[0] <<endl;

	return 0;
}
double apple::getWeight(){
	return weight; 
}
	

