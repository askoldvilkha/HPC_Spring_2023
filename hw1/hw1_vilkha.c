#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*this function is the function we want to integrate
it takes as input:
x - the variable*/
double func(double x) 
{
    return sqrt(1-x*x);
}

/*this function is filling our array 
it "tabulates" the x and f(x) for the following integration
it takes as input:
f_array - array to be filled
n - number of elements
a - lower limit of the integral
dx - "step" of integration
*/
void array_filler(double *f_array, int n, double a, double dx)
{
    double x_i;

    for (int i=0; i<n+1; i++)
    {
        x_i=a+i*dx;
        f_array[i]=func(x_i);
    }
}

/*this function is implementing the simpson's rule
it takes as input:
f_array - array of tabulated! function values
n - number of steps
dx - step of integration
*/
double simpson(double *f_array, int n, double dx)
{

    double sum=0;
    sum=f_array[0]-f_array[n];
    // YZ: If you split into two loops, you can remove all branches inside
    // the loop
    for (int i=1; i<n; i=i+2)
    {
        sum=sum+4*f_array[i]+2*f_array[i+1];
    }
    sum=sum*dx/3;
    
    return sum;
}

int main(int argc, char **argv)
{
    /*here we check if the input is valid and ask to enter correct value if it's not*/
    // YZ: You used agv[1] before ensuring that it exists
    if (argc!=2)
    {
        printf("Number of arguments is not valid! Enter 1 positive and even number next time \n");
        exit(-1);
    }
    int n=atoi(argv[1]);
    // YZ: while (n==n) looks strange. Since you "want" an infinite
    // loop, use while (1). Also, I suggest just reporting an error
    // and quiting if a user supplies an invalid inpu
    /*this "infinite loop" stops whenever user enters a valid input*/
    while (1)
    {
        char input[10];
        if(n < 2)
        {
        printf("You either entered a negative number or an infinitely large number. \n");
        printf("Both of these are not valid and will produce bugs. Try a different positive and even number. \n");
        printf("Note! The program will read only the first number! \n");
        printf("For your reference - at 1000 steps error is of -1.654232e-14. You don't need so many steps! \n");
        scanf("%s", input);
        }
        else if(n % 2 != 0)
        {
        printf("That's not an even integer, enter a positive even number \n");
        printf("Note! The program will read only the first number! \n");
        scanf("%s", input);
        }
        else if(n>1024*1024)
        {
        printf("That's too big, enter a smaller number unless you want to get a bug \n");
        printf("Note! The program will read only the first number! \n");
        scanf("%s", input);
        }
        else 
        {
            break;
        }
        n=atoi(input);
    }
    /*allocating memory*/
    double *f_array = malloc((n+1)*sizeof(double));
    /*defining limits of the integral and step of integration*/
    double a=0;
    double b=1/sqrt(2);
    double dx=(b-a)/n;
    /*filling the array ie tabulating the function*/
    array_filler(f_array, n, a, dx);

    double sum_n=simpson(f_array, n, dx);
    /*calculating the real value of the integral to compare and find the error*/
    double res;
    double pi=3.14159265358979323846;
    res=pi/8+.25*sin(pi/2);
    double err=sum_n-res;
    /*producing output*/
    printf("%e %e %e \n", dx, sum_n, err); 

    return 0;

}
