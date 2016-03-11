#include <FPT.h>
#include <D3d_matrix.h>
#include <xwd_tools.h>

// tetrahedron model
char textureName[100];

int sphere(double u, double v, double p[3]){
  p[0] = sqrt(1-(v*v))*cos(u);
  p[1] = v;
  p[2] = sqrt(1-(v*v))*sin(u);
  return 1;
}

int hyperboloid(double u, double v, double p[3]){
  p[0]=sqrt(1+(v*v))*cos(u);
  p[1]=v;
  p[2]=sqrt(1+(v*v))*sin(u);
  return 1;
}

int crossProduct(double a[3], double b[3], double normalVector[3]){
  normalVector[0] = (a[1]*b[2])-(a[2]*b[1]);
  normalVector[1] = (a[2]*b[0])-(a[0]*b[2]);
  normalVector[2] = (a[0]*b[1])-(a[1]*b[0]);
  return 1;
}

double z_Buffer[600][600];

// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;



int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}

//========================================================================

void makeManMatrix(double manMatrix[4][4], double invManMatrix[4][4], double sx, double sy, double sz, double rx, double ry, double rz, double tx, double ty, double tz){
  int num=0;
  int tlist[9];
  double plist[9];

  tlist[num] = SX ; plist[num] =  sx ; num++ ;
  tlist[num] = SY ; plist[num] =  sy ; num++ ;
  tlist[num] = SZ ; plist[num] =  sz ; num++ ;
  tlist[num] = RX ; plist[num] =  rx ; num++ ;
  tlist[num] = RY ; plist[num] =  ry ; num++ ;
  tlist[num] = RZ ; plist[num] =  rz ; num++ ;
  tlist[num] = TX ; plist[num] =  tx ; num++ ;
  tlist[num] = TY ; plist[num] =  ty ; num++ ;
  tlist[num] = TZ ; plist[num] =  tz ; num++ ;


  D3d_make_movement_sequence_matrix(manMatrix, invManMatrix, num, tlist, plist);
}

//==============================================================================



void drawobject(int textureFlag, double eye[3], double halfAngle, double light[3], double view[4][4], double inherentRGB[3], double manMatrix[4][4], int (*func)(double u, double v, double p[3])){

  double p[3] , q[3], r[3], n[3], xBar, yBar, tanHalf, tempP[3];
  tanHalf = tan(halfAngle*(M_PI/180));
  double actualRGB[3], tempX, tempY;
  double i, j, e, width, height;

  int textureMap, d[2];
  if( textureFlag == 1){
    textureMap = init_xwd_map_from_file(textureName);
    if(textureMap == -1){ printf("failure\n"); }
    e = get_xwd_map_dimensions(textureMap, d);
    if(e == -1){ printf("failure2\n");}
    width = d[0]; height = d[1];
    printf("width: %lf, height: %lf\n", width, height);
  }
  for(i=0; i<2*M_PI; i+=.001){
    for(j=-1; j<1; j+=.001){
        func(i, j, p);
        D3d_mat_mult_pt(p, manMatrix, p);
        tempP[0] = p[0];
        tempP[1] = p[1];
        tempP[2] = p[2];
        D3d_mat_mult_pt(tempP, view, tempP);
        if(tempP[2] < 0) continue;
        else if(fabs(tempP[1]/tempP[2]) > tanHalf) continue;
        else if(fabs(tempP[0]/tempP[2]) > tanHalf) continue;
        light_in_eye_space[0] = light[0] ;
        light_in_eye_space[1] = light[1] ;
        light_in_eye_space[2] = light[2] ;
        //light model stuff
        func(i+.01, j, q);
        func(i, j+.01, r);
        D3d_mat_mult_pt(q, manMatrix, q);
        D3d_mat_mult_pt(r, manMatrix, r);
        q[0] = q[0]-p[0]; q[1] = q[1]-p[1]; q[2] = q[2]-p[2];
        r[0] = r[0]-p[0]; r[1] = r[1]-p[1]; r[2] = r[2]-p[2];
        crossProduct(r, q, n);

        if(textureFlag == 1){
          e = get_xwd_map_color(textureMap, (int)xBar, (int)yBar, inherentRGB);
          if(e==-1){ printf("failure3\n");}
        }

        Light_Model(inherentRGB, eye, p, n, actualRGB);
        G_rgb(actualRGB[0], actualRGB[1], actualRGB[2]);

        D3d_mat_mult_pt(p, view, p);
        xBar = ((300*p[0])/(tanHalf*p[2]))+300;
        yBar = ((300*p[1])/(tanHalf*p[2]))+300;
        //checking the z_Buffer
        if(z_Buffer[(int)xBar][(int)yBar] > p[2]){
          z_Buffer[(int)xBar][(int)yBar] = p[2];
          G_point(xBar, yBar);
        }
      }
    }
}

//====================================================================

void path (int frame_number, double path_xyz[3])
{
  double u,v,r ;
  double x,y,z ;

  u = 5*frame_number*M_PI/180 ;
  v = 0.3*u ;
  r = 2.0 + 1.4*sin(u) ;

  x = r*cos(u)*cos(v) ;
  y = r*sin(u) ;
  z = r*cos(u)*sin(v) ;
 
  path_xyz[0] = x ;
  path_xyz[1] = y ;
  path_xyz[2] = z ;

}

//=================================================================================

int init_scene (int frame_number)
{
  // model variables
  double xcen[4],ycen[4],zcen[4],brad ; // four nodes of tetrahedron
  double ccx,ccy,ccz,ccr ; // location of center of center sphere and radius

  double degrees_of_half_angle = 25;
  double eye[3],coi[3],up[3] ;
  double light_position[3], amb, diff, spow ;
  double theta;
  int k, i, j;
  double V[4][4], Vi[4][4], manMatrix[4][4], invManMatrix[4][4], inherentRGB[3];
  double tempEye[3], tempCoi[3], tempUp[3], Vtemp[4][4], Vitemp[4][4], Len;


  //////////////////////////////////////////////
  //////////////////////////////////////////////
  // build a ball and stick model of a tetrahedron
  //////////////////////////////////////////////
  //////////////////////////////////////////////
  // 3 equally spaced pts around unit circle in the xz-plane 
  // form the base

  for (k = 0 ; k < 3 ; k++) {
    theta = 2*M_PI*k/3 ;
    xcen[k] = cos(theta) ;
    ycen[k] = 0 ;
    zcen[k] = sin(theta) ;
  }

  // you figure where the 4th node of the regular tetrahedron
  xcen[3] = 0 ; ycen[3] = 1 ; zcen[3] = 0 ;

  // also, figure out location of the 5th node of the model
  // which is at the center of mass of the tetrahedron
  for(i=0; i<4; i++){
    ccx = xcen[i];
    ccy = ycen[i]; 
    ccz = zcen[i];
  }

  ccx=ccx/4;
  ccy=ccy/4;
  ccz=ccz/4;

  brad = 0.08 ; // radius of the 4 verts of the tetrahedron
  ccr  = 0.20 ; // the radius of the center node of the model


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  path (frame_number, eye) ;


  coi[0] = ccx ;
  coi[1] = ccy ;
  coi[2] = ccz ;

  path (frame_number + 1, up) ;

  // printf("eye = %lf %lf %lf\n",eye[0],eye[1],eye[2]) ;
  // printf("coi = %lf %lf %lf\n",coi[0],coi[1],coi[2]) ;
  // printf("up  = %lf %lf %lf\n",up[0],up[1],up[2]) ;

  //////////////////////////////////////////////
  //////////////////////////////////////////////
  for(i=0; i<600; i++){
    for(j=0; j<600; j++){
      z_Buffer[i][j]=100000000;
    }
  }

  path (frame_number + 10, light_position) ;
  amb  = 0.2 ;
  diff = 0.5 ;
  spow = 80 ;

  D3d_view(V, Vi, eye, coi, up);

  //center sphere
  inherentRGB[0]=.2; inherentRGB[1]=.4; inherentRGB[2]=.4;
  makeManMatrix(manMatrix, invManMatrix, ccr, ccr, ccr, 0, 0, 0, ccx, ccy, ccz);
  drawobject(1, eye, degrees_of_half_angle, light_position, V, inherentRGB, manMatrix, sphere);
  //verts of tetrahedron
  for(i=0; i<4; i++){
    //printf("xcen: %lf, ycen: %lf, zcen: %lf\n", xcen[i], ycen[i], zcen[i]);
    inherentRGB[0]=.8; inherentRGB[1]=.1; inherentRGB[2]=.1;
    makeManMatrix(manMatrix, invManMatrix, brad, brad, brad, 0, 0, 0, xcen[i], ycen[i], zcen[i]);
    drawobject(0, eye, degrees_of_half_angle, light_position, V, inherentRGB, manMatrix, sphere);
  }
  //translate hyperboloid up, then rotate to lay flat on z axis
  //use view inverse, with eye as starting sphere, and coi as ending sphere
  //hyperboloids
  inherentRGB[0]=.1; inherentRGB[1]=.8; inherentRGB[2]=.1;
  //sx, sy, sz, rx, ry, rz, tx, ty, tz
  for(i=0; i<4; i++){
    for(k=i; k<4; k++){
      if(i==k) continue;
        tempEye[0] = xcen[i];
        tempEye[1] = ycen[i];
        tempEye[2] = zcen[i];

        tempCoi[0] = xcen[k];
        tempCoi[1] = ycen[k];
        tempCoi[2] = zcen[k];

        tempUp[0] = xcen[i];
        tempUp[1] = ycen[i]+1;
        tempUp[2] = zcen[i];
        Len = sqrt(pow(xcen[k]-xcen[i],2) + pow(ycen[k]-ycen[i],2) + pow(zcen[k]-zcen[i],2)) ;
        D3d_view(Vtemp, Vitemp, tempEye, tempCoi, tempUp);
        makeManMatrix(manMatrix, invManMatrix, .03, Len/2, .03,   90, 0, 0,   0, 0, Len/2);
        D3d_mat_mult(manMatrix, Vitemp, manMatrix);
        drawobject(0, eye, degrees_of_half_angle, light_position, V, inherentRGB, manMatrix, hyperboloid);
    }
  }
//hyperboloids that connect center to verts
  tempEye[0] = ccx;
  tempEye[1] = ccy;
  tempEye[2] = ccz;

  tempUp[0] = ccx+.5;
  tempUp[1] = ccy+.3;
  tempUp[2] = ccz+.8;
  inherentRGB[0] = .1; inherentRGB[1] = .1; inherentRGB[2] = .8;
  for(i=0; i<4; i++){
    tempCoi[0] = xcen[i];
    tempCoi[1] = ycen[i];
    tempCoi[2] = zcen[i];
    Len = sqrt(pow(ccx-xcen[i],2) + pow(ccy-ycen[i],2) + pow(ccz-zcen[i],2)) ;
    D3d_view(Vtemp, Vitemp, tempEye, tempCoi, tempUp);
    makeManMatrix(manMatrix, invManMatrix, .02, Len/2, .02,   90, 0, 0,   0, 0, Len/2);
    D3d_mat_mult(manMatrix, Vitemp, manMatrix);
    drawobject(0, eye, degrees_of_half_angle, light_position, V, inherentRGB, manMatrix, hyperboloid);
  }
}

int main(){
  int x;
  G_init_graphics(600,600);
  int i;
  char prefix[100], filename[100];
  strncpy(prefix, "tetrahedron", 100);
  scanf("%s", textureName);

  //tetrahedron0000
  //14
  for(i=0; i < 73; i++){
    printf("%d\n", i);
    init_scene(i);
    sprintf(filename, "%s%04d", prefix, i);
    filename[15] = '.';
    filename[16] = 'x';
    filename[17] = 'w';
    filename[18] = 'd';
    G_save_image_to_file(filename);
    G_wait_key();
    G_rgb(1,1,1);
    G_clear();
  }
//main->init scene->draw
}


