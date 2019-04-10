#include <iostream>
#include "Matrix.h"
using namespace std;

int main()
{
	Matrix m(3,3);
	cin >> m;
	cout << m.inverse() << endl;
	cout << m*m.inverse();

	/*try 
	{
		cout<<m*n;
	}
	catch(...)
	{
		cout << "error";
	}*/
	system("pause");
	return 0;
}