#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415927
#endif

fwd_kin(theta, x)
double *theta;
double x[3];
{
    double l[4];
    l[0]=0.25;
    l[1]=0.2;
    l[2]=0.2;
    l[3]=0.15;
    
    double d[5];
    d[0]=0; //extra created to match d1, d2, d3, d4 as given in assignment
    d[1]=-0.04;
    d[2]=0.04;
    d[3]=-0.04;
    d[4]=-0.04;
    
    double Dx_l1[4][4];
    double Dx_l2[4][4];
    double Dx_l3[4][4];
    
    double Dy_d1[4][4];
    double Dy_d2[4][4];
    double Dy_d3[4][4];
    
    double Dz_l0[4][4];
    double Dz_d4[4][4];

    double Rx_theta0[4][4];
    double Rx_theta1[4][4];
    double Rx_theta2[4][4];
    double Rx_theta3[4][4];

    double Ry_theta0[4][4];
    double Ry_theta1[4][4];
    double Ry_theta2[4][4];
    double Ry_theta3[4][4];

    double Rz_theta0[4][4];
    double Rz_theta1[4][4];
    double Rz_theta2[4][4];
    double Rz_theta3[4][4];

    double sin(double redian);
    double cos(double rad);
    
    double Transformation_Matrix[4][4];

    // fill translational matrices
    for (int row=0; row<4; row++){
        for (int col=0; col<4; col++){
            if (row==col){
                Dx_l1[row][col]=1;
                Dx_l2[row][col]=1;
                Dx_l3[row][col]=1;
                
                Dy_d1[row][col]=1;
                Dy_d2[row][col]=1;
                Dy_d3[row][col]=1;
                
                Dz_l0[row][col]=1;
                Dz_d4[row][col]=1;
            }
            else if(row==0 && col==3){
                 Dx_l1[row][col]=l[1];
                 Dx_l2[row][col]=l[2];
                 Dx_l3[row][col]=l[3];
                 
                 Dy_d1[row][col]=0;
                 Dy_d2[row][col]=0;
                 Dy_d3[row][col]=0;
                 
                 Dz_l0[row][col]=0;
                 Dz_d4[row][col]=0;
            }
            else if(row==1 && col==3){
                 Dx_l1[row][col]=0;
                 Dx_l2[row][col]=0;
                 Dx_l3[row][col]=0;
                 
                 Dy_d1[row][col]=d[1];
                 Dy_d2[row][col]=d[2];
                 Dy_d3[row][col]=d[3];
                 
                 Dz_l0[row][col]=0;
                 Dz_d4[row][col]=0;
            }
            else if(row==2 && col==3){
                 Dx_l1[row][col]=0;
                 Dx_l2[row][col]=0;
                 Dx_l3[row][col]=0;
                 
                 Dy_d1[row][col]=0;
                 Dy_d2[row][col]=0;
                 Dy_d3[row][col]=0;
                 
                 Dz_l0[row][col]=l[0];
                 Dz_d4[row][col]=d[4];
            }
            else {
                 Dx_l1[row][col]=0;
                 Dx_l2[row][col]=0;
                 Dx_l3[row][col]=0;
                 
                 Dy_d1[row][col]=0;
                 Dy_d2[row][col]=0;
                 Dy_d3[row][col]=0;
                 
                 Dz_l0[row][col]=0;
                 Dz_d4[row][col]=0;
            }
        }
    }

    // fill rotational matrices
    // -------------------------------------fill Rx matrices------------------------------------
    // fill Rx_theta0
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if (r==0 && c==0){
                Rx_theta0[r][c]=1;
            }
            else if(r==1 && c==1){
                Rx_theta0[r][c]=cos(theta[0]);
            }
            else if(r==2 && c==1){
                Rx_theta0[r][c]=sin(theta[0]);
            }
            else if(r==1 && c==2){
                Rx_theta0[r][c]= -sin(theta[0]);
            }
            else if(r==2 && c==2){
                Rx_theta0[r][c]= cos(theta[0]);
            }
            else if(r==3 && c==3){
                Rx_theta0[r][c]= 1;
            }
            else{
                Rx_theta0[r][c]=0;
            }
        }
    }
    
    // fill Rx_theta1
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if (r==0 && c==0){
                Rx_theta1[r][c]=1;
            }
            else if(r==1 && c==1){
                Rx_theta1[r][c]=cos(theta[1]);
            }
            else if(r==2 && c==1){
                Rx_theta1[r][c]=sin(theta[1]);
            }
            else if(r==1 && c==2){
                Rx_theta1[r][c]= -sin(theta[1]);
            }
            else if(r==2 && c==2){
                Rx_theta1[r][c]= cos(theta[1]);
            }
            else if(r==3 && c==3){
                Rx_theta1[r][c]= 1;
            }
            else{
                Rx_theta1[r][c]=0;
            }
        }
    }
    
    // fill Rx_theta2
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if (r==0 && c==0){
                Rx_theta2[r][c]=1;
            }
            else if(r==1 && c==1){
                Rx_theta2[r][c]=cos(theta[2]);
            }
            else if(r==2 && c==1){
                Rx_theta2[r][c]=sin(theta[2]);
            }
            else if(r==1 && c==2){
                Rx_theta2[r][c]= -sin(theta[2]);
            }
            else if(r==2 && c==2){
                Rx_theta2[r][c]= cos(theta[2]);
            }
            else if(r==3 && c==3){
                Rx_theta2[r][c]= 1;
            }
            else{
                Rx_theta2[r][c]=0;
            }
        }
    }
    
    // fill Rx_theta3
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if (r==0 && c==0){
                Rx_theta3[r][c]=1;
            }
            else if(r==1 && c==1){
                Rx_theta3[r][c]=cos(theta[3]);
            }
            else if(r==2 && c==1){
                Rx_theta3[r][c]=sin(theta[3]);
            }
            else if(r==1 && c==2){
                Rx_theta3[r][c]= -sin(theta[3]);
            }
            else if(r==2 && c==2){
                Rx_theta3[r][c]= cos(theta[3]);
            }
            else if(r==3 && c==3){
                Rx_theta3[r][c]= 1;
            }
            else{
                Rx_theta3[r][c]=0;
            }
        }
    }
    
    // -------------------------------------fill Ry matrices------------------------------------
    // fill Ry_theta0
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Ry_theta0[r][c]=cos(theta[0]);
            }
            else if(r==2 && c==0){
                Ry_theta0[r][c]=-sin(theta[0]);
            }
            else if(r==1 && c==1){
                Ry_theta0[r][c]=1;
            }
            else if(r==0 && c==2){
                Ry_theta0[r][c]=sin(theta[0]);
            }
            else if(r==2 && c==2){
                Ry_theta0[r][c]=cos(theta[0]);
            }
            else if(r==3 && c==3){
                Ry_theta0[r][c]=1;
            }
            else{
                Ry_theta0[r][c]=0;
            }
        }
    }
    
    // fill Ry_theta1
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Ry_theta1[r][c]=cos(theta[1]);
            }
            else if(r==2 && c==0){
                Ry_theta1[r][c]=-sin(theta[1]);
            }
            else if(r==1 && c==1){
                Ry_theta1[r][c]=1;
            }
            else if(r==0 && c==2){
                Ry_theta1[r][c]=sin(theta[1]);
            }
            else if(r==2 && c==2){
                Ry_theta1[r][c]=cos(theta[1]);
            }
            else if(r==3 && c==3){
                Ry_theta1[r][c]=1;
            }
            else{
                Ry_theta1[r][c]=0;
            }
        }
    }
    
    // fill Ry_theta2
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Ry_theta2[r][c]=cos(theta[2]);
            }
            else if(r==2 && c==0){
                Ry_theta2[r][c]=-sin(theta[2]);
            }
            else if(r==1 && c==1){
                Ry_theta2[r][c]=1;
            }
            else if(r==0 && c==2){
                Ry_theta2[r][c]=sin(theta[2]);
            }
            else if(r==2 && c==2){
                Ry_theta2[r][c]=cos(theta[2]);
            }
            else if(r==3 && c==3){
                Ry_theta2[r][c]=1;
            }
            else{
                Ry_theta2[r][c]=0;
            }
        }
    }
    
    // fill Ry_theta3
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Ry_theta3[r][c]=cos(theta[3]);
            }
            else if(r==2 && c==0){
                Ry_theta3[r][c]=-sin(theta[3]);
            }
            else if(r==1 && c==1){
                Ry_theta3[r][c]=1;
            }
            else if(r==0 && c==2){
                Ry_theta3[r][c]=sin(theta[3]);
            }
            else if(r==2 && c==2){
                Ry_theta3[r][c]=cos(theta[3]);
            }
            else if(r==3 && c==3){
                Ry_theta3[r][c]=1;
            }
            else{
                Ry_theta3[r][c]=0;
            }
        }
    }
    
    // -------------------------------------fill Rz matrices------------------------------------
    // fill Rz_theta0
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Rz_theta0[r][c]=cos(theta[0]);
            }
            else if(r==1 && c==0){
                Rz_theta0[r][c]=sin(theta[0]);
            }
            else if (r==0 && c==1){
                Rz_theta0[r][c]=-sin(theta[0]);
            }
            else if (r==1 && c==1){
                Rz_theta0[r][c]=cos(theta[0]);
            }
            else if (r==2 && c==2){
                Rz_theta0[r][c]=1;
            }
            else if (r==3 && c==3){
                Rz_theta0[r][c]=1;
            }
            else{
                Rz_theta0[r][c]=0;
            }
        }
    }
    
    // fill Rz_theta1
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Rz_theta1[r][c]=cos(theta[1]);
            }
            else if(r==1 && c==0){
                Rz_theta1[r][c]=sin(theta[1]);
            }
            else if (r==0 && c==1){
                Rz_theta1[r][c]=-sin(theta[1]);
            }
            else if (r==1 && c==1){
                Rz_theta1[r][c]=cos(theta[1]);
            }
            else if (r==2 && c==2){
                Rz_theta1[r][c]=1;
            }
            else if (r==3 && c==3){
                Rz_theta1[r][c]=1;
            }
            else{
                Rz_theta1[r][c]=0;
            }
        }
    }
    
    // fill Rz_theta2
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Rz_theta2[r][c]=cos(theta[2]);
            }
            else if(r==1 && c==0){
                Rz_theta2[r][c]=sin(theta[2]);
            }
            else if (r==0 && c==1){
                Rz_theta2[r][c]=-sin(theta[2]);
            }
            else if (r==1 && c==1){
                Rz_theta2[r][c]=cos(theta[2]);
            }
            else if (r==2 && c==2){
                Rz_theta2[r][c]=1;
            }
            else if (r==3 && c==3){
                Rz_theta2[r][c]=1;
            }
            else{
                Rz_theta2[r][c]=0;
            }
        }
    }
    
    // fill Rz_theta3
    for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
            if(r==0 && c==0){
                Rz_theta3[r][c]=cos(theta[3]);
            }
            else if(r==1 && c==0){
                Rz_theta3[r][c]=sin(theta[3]);
            }
            else if (r==0 && c==1){
                Rz_theta3[r][c]=-sin(theta[3]);
            }
            else if (r==1 && c==1){
                Rz_theta3[r][c]=cos(theta[3]);
            }
            else if (r==2 && c==2){
                Rz_theta3[r][c]=1;
            }
            else if (r==3 && c==3){
                Rz_theta3[r][c]=1;
            }
            else{
                Rz_theta3[r][c]=0;
            }
        }
    }
    
    // Matrices multiplications
    //Transformation_Matrix[4][4] = Dz_l0*Rz_theta0*Ry_theta1*Dy_d1*Dx_l1*Ry_theta2*Dy_d2*Dx_l2*Ry_theta3*Dy_d3*Dz_d4*Dx_l3; //reverse it
        
    // matrices multiplication (A * B)
    int n=4; // 4*4 matrices
    double temp = 0;

    //Dz_d4*Dx_l3
    double T10[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Dz_d4[rowA][colA] * Dx_l3[colA][colB];
            }
            T10[rowA][colB] = temp;
            temp=0;
        }
    }

    //Dy_d3*T10
    double T9[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Dy_d3[rowA][colA] * T10[colA][colB];
            }
            T9[rowA][colB] = temp;
            temp=0;
        }
    }

    //Ry_theta3*T9
    double T8[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Ry_theta3[rowA][colA] * T9[colA][colB];
            }
            T8[rowA][colB] = temp;
            temp=0;
        }
    }

    //Dx_l2*T8
    double T7[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Dx_l2[rowA][colA] * T8[colA][colB];
            }
            T7[rowA][colB] = temp;
            temp=0;
        }
    }

    //Dy_d2*T7
    double T6[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Dy_d2[rowA][colA] * T7[colA][colB];
            }
            T6[rowA][colB] = temp;
            temp=0;
        }
    }

    //Ry_theta2*T6
    double T5[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Ry_theta2[rowA][colA] * T6[colA][colB];
            }
            T5[rowA][colB] = temp;
            temp=0;
        }
    }

    //Dx_l1*T5
    double T4[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Dx_l1[rowA][colA] * T5[colA][colB];
            }
            T4[rowA][colB] = temp;
            temp=0;
        }
    }

    //Dy_d1*T4
    double T3[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Dy_d1[rowA][colA] * T4[colA][colB];
            }
            T3[rowA][colB] = temp;
            temp=0;
        }
    }

    //Ry_theta1*T3
    double T2[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Ry_theta1[rowA][colA] * T3[colA][colB];
            }
            T2[rowA][colB] = temp;
            temp=0;
        }
    }

    //Rz_theta0*T2
    double T1[4][4];
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Rz_theta0[rowA][colA] * T2[colA][colB];
            }
            T1[rowA][colB] = temp;
            temp=0;
        }
    }

    //Dz_l0*T1
    for(int rowA=0; rowA<n; rowA++){
        for(int colB=0; colB<n; colB++){
            for(int colA=0; colA<n; colA++){
                    temp = temp + Dz_l0[rowA][colA] * T1[colA][colB];
            }
            Transformation_Matrix[rowA][colB] = temp;
            temp=0;
        }
    }

    //assign positions of tool frame from transformation matrix obtained
    x[0] = Transformation_Matrix[0][3];
    x[1] = Transformation_Matrix[1][3];
    x[2] = Transformation_Matrix[2][3];

    return *x;
}

inv_kin(x, theta)
double *x;
double theta[6];
{
    double l[4];
    l[0]=0.25;
    l[1]=0.2;
    l[2]=0.2;
    l[3]=0.15;

    double d[5];
    d[0]=0; //extra created to match d1, d2, d3, d4 as given in assignment
    d[1]=-0.04;
    d[2]=0.04;
    d[3]=-0.04;
    d[4]=-0.04;

    double acos(double x);
    double atan(double x);
    double sqrt(double arg);

    double temp =0;
    double temp1 =0;
    double temparg =0;
    
    double alpha;
    double gamma;
    double plustheta0;

    temparg = d[3]/(l[1]+l[2]);
    plustheta0 = atan(temparg);
    temparg = 0;

    double x1[3];// x', y', z'
    x1[0]=x[0]+d[4];
    x1[1]=x[1]; 
    x1[2]=x[2]+l[3];

    double x2[3];// x'', y'', z''
    x2[0]=sqrt(x1[0]*x1[0]+x1[1]*x1[1]);
    x2[1]=x1[1];
    x2[2]=x1[2]-l[0];

    theta[0] = atan(x1[1]/x1[0]) - plustheta0;
    temp = atan(x1[1]/x1[0]) + plustheta0;
    temp = 0;

    // theta2 = (+-)acos((x''^2 +y''^2)-l1^2-l2^2)/2*l1*l2)
    theta[2] = acos((x2[0]*x2[0]+x2[1]*x2[1]-l[1]*l[1]-l[2]*l[2])/(2*l[1]*l[2]));

    // alpha = acos((l1^2+x''^2+y''^2-l2^2)/(2*l1*sqrt(x''^2+y''^2))
    alpha = acos((l[1]*l[1]+x2[0]*x2[0]+x2[1]*x2[1]-l[2]*l[2])/(2*l[1]*sqrt(x2[0]*x2[0]+x2[1]*x2[1])));

    temparg = x2[1]/x2[0];
    gamma = atan(temparg); //gamma = aTan(y''/x'')
    temparg = 0;
    
    temp = gamma + alpha;
    temp1 = gamma - alpha;
    theta[1] = temp;
    temp = 0;
    temp1 = 0;
  
    //from triangles formed by substructure of links
    theta[3] = M_PI + theta[2] + alpha - (M_PI/2) + gamma; // theta3  = theta'(i.e Pi + theta2) + alpha - Pi/2 + gamma

    theta[4] = 0;
    
    return theta;
}