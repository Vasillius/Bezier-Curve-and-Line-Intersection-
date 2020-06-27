/*this code is heavily based on code posted on https://www.particleincell.com/2013/cubic-line-intersection/ , please chech it for for explanation
also check: 
https://www.youtube.com/watch?v=TeXajQ62yZ8 
https://www.youtube.com/watch?v=pnYccz1Ha34
https://www.youtube.com/watch?v=2HvH9cmHbG4&t=1307s
for a better understanding 
*/


int x1, x2; // xs of the endpoints 
int y1, y2; // ys of the end point 
float c1x, c1y, c2x, c2y; //coordinates of the control points 
PVector l1, l2; // end-points of the line 
PVector[] bp; // cubic bezier's four points 
double t; // t value of the bezier curve between (0 - 1)

void setup() {

  size(500, 500);


  //setting up the bezier's points' positions 
  x1 = width / 4;
  y1  = height / 8;

  x2 = width - (width / 4);
  y2  = height / 8;

  c1x = x1 + 10;
  c1y= height/2;

  c2x = x2 - 10;
  c2y = height/2;

  l1 = new PVector (width/2, 0);
  
  bp = new PVector[4]; // stores the four points as vectors 
}



void draw() {


  background(255);
  noFill();
  stroke(0);
  strokeWeight(1);
  bezier(x1, y1, c1x, c1y, c2x, c2y, x2, y2); //drawing the bzier curve 



  bp[0] = new PVector(x1, y1);
  bp[1] = new PVector(c1x, c1y);
  bp[2] = new PVector(c2x, c2y);
  bp[3] = new PVector(x2, y2);


  l2 = new PVector (mouseX, mouseY);// line's end-point is at mouse's position 

  strokeWeight(1);
  line(l1.x, l1.y, l2.x, l2.y);

  computeIntersections(bp, l1, l2);
}





/*computes intersection between a cubic spline and a line segment*/
void computeIntersections(PVector[] p, PVector l1, PVector l2) // takes the 4 points of the bezier curve plus the the two endpoints of the line 
{
  float[] X = new float[2]; // this will store the coordinates of the intersection point 

  float A = l2.y - l1.y;      //A=y2-y1
  float B = l1.x - l2.x;      //B=x1-x2
  float C = l1.x * (l1.y - l2.y) + l1.y * (l2.x - l1.x);//C=x1*(y1-y2)+y1*(x2-x1)

  float [] bx = bezierCoeffs(p[0].x, p[1].x, p[2].x, p[3].x);
  float [] by = bezierCoeffs(p[0].y, p[1].y, p[2].y, p[3].y);

  float[] P = new float[4];
  P[0] = A*bx[0]+B*by[0];    /*t^3*/
  P[1] = A*bx[1]+B*by[1];    /*t^2*/
  P[2] = A*bx[2]+B*by[2];    /*t*/
  P[3] = A*bx[3]+B*by[3] + C;  /*1*/

  float[] r = cubicRoots(P);

  /*verify the roots are in bounds of the linear segment*/
  for (int i = 0; i < 3; i++)
  {
    t = r[i];


    X[0] = (float)(bx[0]*t*t*t+bx[1]*t*t+bx[2]*t+bx[3]);
    X[1] = (float)(by[0]*t*t*t+by[1]*t*t+by[2]*t+by[3]);            

    /*above is intersection point assuming infinitely long line segment,
     make sure we are also in bounds of the line*/
    float s;
    if ((l2.x - l1.x)!=0) {          /*if not vertical line*/
      s = (X[0] - l1.x)/(l2.x - l1.x);
    } else {
      s = (X[1] - l1.y)/(l2.y - l1.y);
    }

    /*in bounds?*/
    if (t<0 || t>1.0 || s<0 || s>1.0)
    {
      X[0]=-100;  /*move off screen*/
      X[1]=-100;
    }


    if (X[0] != -100 && X[1] != -100) {
      fill(255, 0, 0);
      noStroke();
      ellipse((float)X[0], (float)X[1], 10, 10); // drawing an ellipse at the intersection point 
    }





    if (t!= -1) {
      
      // get the coordinates at a specific t value across the curve 
      float x = bezierPoint(x1, c1x, c2x, x2, (float)t);
      float y = bezierPoint(y1, c1y, c2y, y2, (float)t);

     // get the tangent point at a specific t value across the curve 
      float tx = bezierTangent(x1, c1x, c2x, x2, (float)t);
      float ty = bezierTangent(y1, c1y, c2y, y2, (float)t);


      float a = atan2(ty, tx); // get the angle of the tangent line 
      a += PI;
      stroke(255, 102, 0);
      //draw two lines going in opposite directions from the intersection point 
      line(x, y, cos(a)*100 + x, sin(a)*100 + y);
      line(x, y, cos(a)*-100 + x, sin(a)*-100 + y);

      // get the angle of the normal line and draw it 
      float b = atan2(ty, tx);
      b -= HALF_PI;
      line(x, y, cos(b)*200 + x, sin(b)*200 + y);

      
      PVector v1 = new PVector(x - (cos(b)*200 + x), y - (sin(b)*200 + y));//convert the coordinates of normal line into vector form 
      PVector v2 = new PVector(l2.x - l1.x, l2.y - l1.y); //convert the coordinates of intersecting line into vector form 


      float angle = angle(v1, v2); //calculate the angle between the normal line and intersection line 
      float anglef = degrees(angle);
      if (anglef < 180) {
        println((anglef));
      } else {
        println(-1 * (360 - anglef));
      }
    }
  }
}





float[] bezierCoeffs(float P0, float P1, float P2, float P3)
{
  float[] Z = new float[4];
  Z[0] = -P0 + 3*P1 + -3*P2 + P3; 
  Z[1] = 3*P0 - 6*P1 + 3*P2;
  Z[2] = -3*P0 + 3*P1;
  Z[3] = P0;
  return Z;
}



float[] cubicRoots(float[] P)
{
  float a = P[0];
  float b = P[1];
  float c = P[2];
  float d = P[3];

  float A = b / a;
  float B = c / a;
  float C = d / a;

  float Q, R, D, S, T, Im;

  Q = (float)((3*B - Math.pow(A, 2))/9);
  R = (float)((9*A*B - 27*C - 2*Math.pow(A, 3))/54);
  D = (float)(Math.pow(Q, 3) + Math.pow(R, 2));    // polynomial discriminant

  float[] t = new float[3];

  if (D >= 0)                                 // complex or duplicate roots
  {
    S = sgn((R + Math.sqrt(D))*Math.pow(Math.abs(R + Math.sqrt(D)), (1/3)));
    T = sgn((R - Math.sqrt(D))*Math.pow(Math.abs(R - Math.sqrt(D)), (1/3)));

    t[0] = -A/3 + (S + T);                    // real root
    t[1] = -A/3 - (S + T)/2;                  // real part of complex root
    t[2] = -A/3 - (S + T)/2;                  // real part of complex root
    Im = (float)(Math.abs(Math.sqrt(3)*(S - T)/2));    // complex part of root pair   

    /*discard complex roots*/
    if (Im!=0)
    {
      t[1]=-1;
      t[2]=-1;
    }
  } else                                          // distinct real roots
  {
    float th = (float)(Math.acos(R/Math.sqrt(-Math.pow(Q, 3))));

    t[0] = (float)(2*Math.sqrt(-Q)*Math.cos(th/3) - A/3);
    t[1] = (float)(2*Math.sqrt(-Q)*Math.cos((th + 2*Math.PI)/3) - A/3);
    t[2] = (float)(2*Math.sqrt(-Q)*Math.cos((th + 4*Math.PI)/3) - A/3);
    Im = 0.0;
  }

  /*discard out of spec roots*/
  for (int i = 0; i < 3; i++) 
    if (t[i]<0 || t[i]>1.0) t[i]=-1;

  /*sort but place -1 at the end*/
  t = sortSpecial(t);

  //console.log(t[0]+" "+t[1]+" "+t[2]);
  return t;
}




int sgn( double x )
{
  if (x < 0.0) { 
    return -1;
  } else {
    return 1;
  }
}




float [] sortSpecial(float [] a)
{
  boolean flip;
  float temp;

  do {
    flip = false;
    for (int i = 0; i < a.length - 1; i++)
    {
      if ((a[i+1]>=0 && a[i]>a[i+1]) ||
        (a[i]<0 && a[i+1]>=0))
      {
        flip = true;
        temp = a[i];
        a[i] = a[i+1];
        a[i+1] = temp;
      }
    }
  } while (flip);
  return a;
}





float angle(PVector v1, PVector v2) {
  float a = atan2(v2.y, v2.x) - atan2(v1.y, v1.x);
  if (a < 0) a += TWO_PI;
  return a;
}
