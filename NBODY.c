#include <stdlib.h>
#include <math.h>
#include <stdio.h>

//void rk4(float f(float, float, float, float, float, float), float v, float x, float h, int n); comente esto porque estaba jodiendo, malparida funcion de mierda te odio.
void rk4();
float accel_t(float x_t, float x_s, float x_l, float x_a, float G, float m);

float accel_t(float x_t, float x_s, float x_l, float x_a, float G, float m)
{
  return ((G * m)*((x_s -x_t)/pow((x_s-x_t),1.5))*((x_l-x_t)/pow((x_l-x_t),1.5))*((x_a-x_t)/pow((x_a-x_t),1.5)));
}

void rk4()
{
	FILE *in;
	// Indice
	int i;
	int j;
	int n = 2; 

	// Nombre del Archivo
	char filename1[30] = "posiciones.txt";

	// Punteros posicion 
	float x_tp;
	float y_tp;
	float z_tp;
	float x_ap;
	float y_ap;
	float z_ap;
	float x_sp;
	float y_sp;
	float z_sp;
	float x_lp;
	float y_lp;
	float z_lp;
	float x_tf;
	float y_tf;
	float z_tf;
	float x_af;
	float y_af;
	float z_af;
	float x_sf;
	float y_sf;
	float z_sf;
	float x_lf;
	float y_lf;
	float z_lf;

	// Punteros  velocidad y vs y as
	float vx_tp;
	float vy_tp;
	float vz_tp;
	float vx_ap;
	float vy_ap;
	float vz_ap;
	float vx_sp;
	float vy_sp;
	float vz_sp;
	float vx_lp;
	float vy_lp;
	float vz_lp;
	float vx_tf;
	float vy_tf;
	float vz_tf;
	float vx_af;
	float vy_af;
	float vz_af;
	float vx_sf;
	float vy_sf;
	float vz_sf;
	float vx_lf;
	float vy_lf;
	float vz_lf;

	// Memory asignationmation
	// x_t = malloc(sizeof(float));
	// y_t = malloc(sizeof(float));
	// z_t = malloc(sizeof(float));
	// x_a = malloc(sizeof(float));
	// y_a = malloc(sizeof(float));
	// z_a = malloc(sizeof(float));
	// x_s = malloc(sizeof(float));
	// y_s = malloc(sizeof(float));
	// z_s = malloc(sizeof(float));
	// x_l = malloc(sizeof(float));
	// y_l = malloc(sizeof(float));
	// z_l = malloc(sizeof(float));
	// vx_t = malloc(sizeof(float));
	// vy_t = malloc(sizeof(float));
	// vz_t = malloc(sizeof(float));
	// vx_a = malloc(sizeof(float));
	// vy_a = malloc(sizeof(float));
	// vz_a = malloc(sizeof(float));
	// vx_s = malloc(sizeof(float));
	// vy_s = malloc(sizeof(float));
	// vz_s = malloc(sizeof(float));
	// vx_l = malloc(sizeof(float));
	// vy_l = malloc(sizeof(float));
	// vz_l = malloc(sizeof(float));

	// aceleraciones y velocidades para calcular vs y xs
	float ax;
	float ay;
	float az;
	float vx;
	float vy;
	float vz;	

	// Constantes Runge-kutta (si, son tantas, para solo llamar la funcion una vez)
	float ktx1,ktx2,ktx3,ktx4,ltx1,ltx2,ltx3,ltx4;
	float ksx1,ksx2,ksx3,ksx4,lsx1,lsx2,lsx3,lsx4;
	float kax1,kax2,kax3,kax4,lax1,lax2,lax3,lax4;
	float klx1,klx2,klx3,klx4,llx1,llx2,llx3,llx4;
	float kty1,kty2,kty3,kty4,lty1,lty2,lty3,lty4;
	float ksy1,ksy2,ksy3,ksy4,lsy1,lsy2,lsy3,lsy4;
	float kay1,kay2,kay3,kay4,lay1,lay2,lay3,lay4;
	float kly1,kly2,kly3,kly4,lly1,lly2,lly3,lly4; 
	float ktz1,ktz2,ktz3,ktz4,ltz1,ltz2,ltz3,ltz4; 
	float ksz1,ksz2,ksz3,ksz4,lsz1,lsz2,lsz3,lsz4;
	float kaz1,kaz2,kaz3,kaz4,laz1,laz2,laz3,laz4;
	float klz1,klz2,klz3,klz4,llz1,llz2,llz3,llz4;

	// Step 
	float h = 0.001;

	// Constante gravitacional
	const float G = 1.0;

	// Masas para los cuatro cuerpos
	const float m_s = 1.0;
	const float m_t = 1*pow(10,-6);
	const float m_l = 1*pow(3.695,-8);
	const float m_a=5*pow(10,-11);

	// C.I. en x y v para los 4 cuerpos en las 3 dimensiones del espacio
    x_tp = 1.0;
    x_sp = 0.0;
    x_ap = 1.0;
    x_lp = 1.0026; 

    y_tp = 0.0;
    y_ap = 0.0963;
    y_sp = 0.0;
    y_lp = 0.0;     

    z_tp = 0.0;
    z_ap = 0.0;
    z_sp = 0.0;
    z_lp = 0.0;

	vx_tp = 0.0;
	vx_sp = 0.0;
	vx_ap = 1570.0;
	vx_lp = 0.0;

	vy_tp = 30.0;
	vy_ap = 0.0;
	vy_sp = 0.0;
	vy_lp = 1.023;

	vz_tp = 0.0;
	vz_ap = 29959.0;
	vz_sp = 0.0;	
	vz_lp = 0.0;
	in = fopen(filename1, "w");
	fprintf(in,"%f \t %f \t %f\n",x_tp,y_tp,z_tp);
	fclose(in);

	in = fopen(filename1, "w");
	for (i = 1; i <= n; i++)
	{
		//Runge-Kutta para x

		// Aqui calculamos k1 para xp y l1 para vxp en los cuatro cuerpos
		ax = accel_t(x_tp,x_sp,x_lp,x_ap,G,m_t);
		vx = vx_tp;
		ktx1 = 0.5 * h * vx;
		ltx1 = 0.5 * h * ax;
		ax = accel_t(x_sp,x_tp,x_lp,x_ap,G,m_s);
		vx = vx_sp;
		ksx1 = 0.5 * h * vx;
		lsx1 = 0.5 * h * ax;
		ax = accel_t(x_ap,x_tp,x_lp,x_sp,G,m_a);
		vx = vx_ap;
		kax1 = 0.5 * h * vx;
		lax1 = 0.5 * h * ax;
		ax = accel_t(x_lp,x_tp,x_sp,x_ap,G,m_l);
		vx = vx_lp;
		klx1 = 0.5 * h * vx;
		llx1 = 0.5 * h * ax;

		// Aqui calculamos k2 para xp y l2 para vxp en los cuatro cuerpos
		ax = accel_t(x_tp+ktx1,x_sp+ksx1,x_lp+klx1,x_ap+kax1,G,m_t);
		vx = vx_tp + ltx1;
		ktx2 = 0.5 * h * vx;
		ltx2 = 0.5 * h * ax;
		ax = accel_t(x_sp+ksx1,x_tp+ktx1,x_lp+klx1,x_ap+kax1,G,m_t);
		vx = vx_sp + lsx1;
		ksx2 = 0.5 * h * vx;
		lsx2 = 0.5 * h * ax;
		ax = accel_t(x_ap+kax1,x_tp+ktx1,x_lp+klx1,x_sp+ksx1,G,m_t);
		vx = vx_ap + lax1;
		kax2 = 0.5 * h * vx;
		lax2 = 0.5 * h * ax;
		ax = accel_t(x_lp+klx1,x_tp+ktx1,x_sp+ksx1,x_ap+kax1,G,m_t);
		vx = vx_lp + llx1;
		klx2 = 0.5 * h * vx;
		llx2 = 0.5 * h * ax;

		// Aqui calculamos k3 para xp y l3 para vxp en los cuatro cuerpos
		ax = accel_t(x_tp+ktx2,x_sp+ksx2,x_lp+klx2,x_ap+kax2,G,m_t);
		vx = vx_tp + ltx2;
		ktx3 = h * vx;
		ltx3 = h * ax;
		ax = accel_t(x_sp+ksx2,x_tp+ktx2,x_lp+klx2,x_ap+kax2,G,m_t);
		vx = vx_sp + lsx2;
		ksx3 = h * vx;
		lsx3 = h * ax;
		ax = accel_t(x_ap+kax2,x_tp+ktx2,x_lp+klx2,x_sp+ksx2,G,m_t);
		vx = vx_ap + lax2;
		kax3 = h * vx;
		lax3 = h * ax;
		ax = accel_t(x_lp+klx2,x_tp+ktx2,x_sp+ksx2,x_ap+kax2,G,m_t);
		vx = vx_lp + llx2;
		klx3 = h * vx;
		llx3 = h * ax;

		// Aqui calculamos k4 para xp y l4 para vxp en los cuatro cuerpos
		ax = accel_t(x_tp+ktx3,x_sp+ksx3,x_lp+klx3,x_ap+kax3,G,m_t);
		vx = vx_tp + ltx3;
		ktx4 = 0.5 * h * vx;
		ltx4 = 0.5 * h * ax;
		ax = accel_t(x_sp+ksx3,x_tp+ktx3,x_lp+klx3,x_ap+kax3,G,m_t);
		vx = vx_sp + lsx3;
		ksx4 = 0.5 * h * vx;
		lsx4 = 0.5 * h * ax;
		ax = accel_t(x_ap+kax3,x_tp+ktx3,x_lp+klx3,x_sp+ksx3,G,m_t);
		vx = vx_ap + lax3;
		kax4 = 0.5 * h * vx;
		lax4 = 0.5 * h * ax;
		ax = accel_t(x_lp+klx3,x_tp+ktx3,x_sp+ksx3,x_ap+kax3,G,m_t);
		vx = vx_lp + llx3;
		klx4 = 0.5 * h * vx;
		llx4 = 0.5 * h * ax;

		//Runge-Kutta para y

		// Aqui calculamos k1 para yp y l1 para vyp en los cuatro cuerpos
		ay = accel_t(y_tp,y_sp,y_lp,y_ap,G,m_t);
		vy = vy_tp;
		kty1 = 0.5 * h * vy;
		lty1 = 0.5 * h * ay;
		ay = accel_t(y_sp,y_tp,y_lp,y_ap,G,m_s);
		vy = vy_sp;
		ksy1 = 0.5 * h * vy;
		lsy1 = 0.5 * h * ay;
		ay = accel_t(y_ap,y_tp,y_lp,y_sp,G,m_a);
		vy = vy_ap;
		kay1 = 0.5 * h * vy;
		lay1 = 0.5 * h * ay;
		ay = accel_t(y_lp,y_tp,y_sp,y_ap,G,m_l);
		vy = vy_lp;
		kly1 = 0.5 * h * vy;
		lly1 = 0.5 * h * ay;

		// Aqui calculamos k2 para yp y l2 para vyp en los cuatro cuerpos
		ay = accel_t(y_tp+kty1,y_sp+ksy1,y_lp+kly1,y_ap+kay1,G,m_t);
		vy = vy_tp + lty1;
		kty2 = 0.5 * h * vy;
		lty2 = 0.5 * h * ay;
		ay = accel_t(y_sp+ksy1,y_tp+kty1,y_lp+kly1,y_ap+kay1,G,m_t);
		vy = vy_sp + lsy1;
		ksy2 = 0.5 * h * vy;
		lsy2 = 0.5 * h * ay;
		ay = accel_t(y_ap+kay1,y_tp+kty1,y_lp+kly1,y_sp+ksy1,G,m_t);
		vy = vy_ap + lay1;
		kay2 = 0.5 * h * vy;
		lay2 = 0.5 * h * ay;
		ay = accel_t(y_lp+kly1,y_tp+kty1,y_sp+ksy1,y_ap+kay1,G,m_t);
		vy = vy_lp + lly1;
		kly2 = 0.5 * h * vy;
		lly2 = 0.5 * h * ay;

		// Aqui calculamos k3 para yp y l3 para vyp en los cuatro cuerpos
		ay = accel_t(y_tp+kty2,y_sp+ksy2,y_lp+kly2,y_ap+kay2,G,m_t);
		vy = vy_tp + lty2;
		kty3 = h * vy;
		lty3 = h * ay;
		ay = accel_t(y_sp+ksy2,y_tp+kty2,y_lp+kly2,y_ap+kay2,G,m_t);
		vy = vy_sp + lsy2;
		ksy3 = h * vy;
		lsy3 = h * ay;
		ay = accel_t(y_ap+kay2,y_tp+kty2,y_lp+kly2,y_sp+ksy2,G,m_t);
		vy = vy_ap + lay2;
		kay3 = h * vy;
		lay3 = h * ay;
		ay = accel_t(y_lp+kly2,y_tp+kty2,y_sp+ksy2,y_ap+kay2,G,m_t);
		vy = vy_lp + lly2;
		kly3 = h * vy;
		lly3 = h * ay;

		// Aqui calculamos k4 para yp y l4 para vyp en los cuatro cuerpos
		ay = accel_t(y_tp+kty3,y_sp+ksy3,y_lp+kly3,y_ap+kay3,G,m_t);
		vy = vy_tp + lty3;
		kty4 = 0.5 * h * vy;
		lty4 = 0.5 * h * ay;
		ay = accel_t(y_sp+ksy3,y_tp+kty3,y_lp+kly3,y_ap+kay3,G,m_t);
		vy = vy_sp + lsy3;
		ksy4 = 0.5 * h * vy;
		lsy4 = 0.5 * h * ay;
		ay = accel_t(y_ap+kay3,y_tp+kty3,y_lp+kly3,y_sp+ksy3,G,m_t);
		vy = vy_ap + lay3;
		kay4 = 0.5 * h * vy;
		lay4 = 0.5 * h * ay;
		ay = accel_t(y_lp+kly3,y_tp+kty3,y_sp+ksy3,y_ap+kay3,G,m_t);
		vy = vy_lp + lly3;
		kly4 = 0.5 * h * vy;
		lly4 = 0.5 * h * ay;

		//Runge Kutta para z

		// Aqui calculamos k1 para zp y l1 para vzp en los cuatro cuerpos
		az = accel_t(z_tp,z_sp,z_lp,z_ap,G,m_t);
		vz = vz_tp;
		ktz1 = 0.5 * h * vz;
		ltz1 = 0.5 * h * az;
		az = accel_t(z_sp,z_tp,z_lp,z_ap,G,m_s);
		vz = vz_sp;
		ksz1 = 0.5 * h * vz;
		lsz1 = 0.5 * h * az;
		az = accel_t(z_ap,z_tp,z_lp,z_sp,G,m_a);
		vz = vz_ap;
		kaz1 = 0.5 * h * vz;
		laz1 = 0.5 * h * az;
		az = accel_t(z_lp,z_tp,z_sp,z_ap,G,m_l);
		vz = vz_lp;
		klz1 = 0.5 * h * vz;
		llz1 = 0.5 * h * az;

		// Aqui calculamos k2 para zp y l2 para vzp en los cuatro cuerpos
		az = accel_t(z_tp+ktz1,z_sp+ksz1,z_lp+klz1,z_ap+kaz1,G,m_t);
		vz = vz_tp + ltz1;
		ktz2 = 0.5 * h * vz;
		ltz2 = 0.5 * h * az;
		az = accel_t(z_sp+ksz1,z_tp+ktz1,z_lp+klz1,z_ap+kaz1,G,m_t);
		vz = vz_sp + lsz1;
		ksz2 = 0.5 * h * vz;
		lsz2 = 0.5 * h * az;
		az = accel_t(z_ap+kaz1,z_tp+ktz1,z_lp+klz1,z_sp+ksz1,G,m_t);
		vz = vz_ap + laz1;
		kaz2 = 0.5 * h * vz;
		laz2 = 0.5 * h * az;
		az = accel_t(z_lp+klz1,z_tp+ktz1,z_sp+ksz1,z_ap+kaz1,G,m_t);
		vz = vz_lp + llz1;
		klz2 = 0.5 * h * vz;
		llz2 = 0.5 * h * az;

		// Aqui calculamos k3 para zp y l3 para vzp en los cuatro cuerpos
		az = accel_t(z_tp+ktz2,z_sp+ksz2,z_lp+klz2,z_ap+kaz2,G,m_t);
		vz = vz_tp + ltz2;
		ktz3 = h * vz;
		ltz3 = h * az;
		az = accel_t(z_sp+ksz2,z_tp+ktz2,z_lp+klz2,z_ap+kaz2,G,m_t);
		vz = vz_sp + lsz2;
		ksz3 = h * vz;
		lsz3 = h * az;
		az = accel_t(z_ap+kaz2,z_tp+ktz2,z_lp+klz2,z_sp+ksz2,G,m_t);
		vz = vz_ap + laz2;
		kaz3 = h * vz;
		laz3 = h * az;
		az = accel_t(z_lp+klz2,z_tp+ktz2,z_sp+ksz2,z_ap+kaz2,G,m_t);
		vz = vz_lp + llz2;
		klz3 = h * vz;
		llz3 = h * az;

		// Aqui calculamos k4 para zp y l4 para vzp en los cuatro cuerpos
		az = accel_t(z_tp+ktz3,z_sp+ksz3,z_lp+klz3,z_ap+kaz3,G,m_t);
		vz = vz_tp + ltz3;
		ktz4 = 0.5 * h * vz;
		ltz4 = 0.5 * h * az;
		az = accel_t(z_sp+ksz3,z_tp+ktz3,z_lp+klz3,z_ap+kaz3,G,m_t);
		vz = vz_sp + lsz3;
		ksz4 = 0.5 * h * vz;
		lsz4 = 0.5 * h * az;
		az = accel_t(z_ap+kaz3,z_tp+ktz3,z_lp+klz3,z_sp+ksz3,G,m_t);
		vz = vz_ap + laz3;
		kaz4 = 0.5 * h * vz;
		laz4 = 0.5 * h * az;
		az = accel_t(z_lp+klz3,z_tp+ktz3,z_sp+ksz3,z_ap+kaz3,G,m_t);
		vz = vz_lp + llz3;
		klz4 = 0.5 * h * vz;
		llz4 = 0.5 * h * az;

		// actualizacion de los vectores con xp hallamos xf 
		/* clarificacion, si ayuda, de alguna manera para nosotros xp es present y xf es future, vivimos en el pasado. */

		x_tf = x_tp + (ktx1 + 2 * (ktx2+ktx3) + ktx4) * 0.16666666666;
		x_sf = x_sp + (ksx1 + 2 * (ksx2+ksx3) + ksx4) * 0.16666666666;
		x_lf = x_lp + (klx1 + 2 * (klx2+klx3) + klx4) * 0.16666666666;
		x_af = x_ap + (kax1 + 2 * (kax2+kax3) + kax4) * 0.16666666666;

		vx_tf = vx_tp + (ltx1 + 2 * (ltx2+ltx3) + ltx4) * 0.16666666666;
		vx_sf = vx_sp + (lsx1 + 2 * (lsx2+lsx3) + lsx4) * 0.16666666666;
		vx_lf = vx_lp + (llx1 + 2 * (llx2+llx3) + llx4) * 0.16666666666;
		vx_af = vx_ap + (lax1 + 2 * (lax2+lax3) + lax4) * 0.16666666666;

		y_tf = y_tp + (kty1 + 2 * (kty2+kty3) + kty4) * 0.16666666666;
		y_sf = y_sp + (ksy1 + 2 * (ksy2+ksy3) + ksy4) * 0.16666666666;
		y_lf = y_lp + (kly1 + 2 * (kly2+kly3) + kly4) * 0.16666666666;
		y_af = y_ap + (kay1 + 2 * (kay2+kay3) + kay4) * 0.16666666666;

		vy_tf = vy_tp + (lty1 + 2 * (lty2+lty3) + lty4) * 0.16666666666;
		vy_sf = vy_sp + (lsy1 + 2 * (lsy2+lsy3) + lsy4) * 0.16666666666;
		vy_lf = vy_lp + (lly1 + 2 * (lly2+lly3) + lly4) * 0.16666666666;
		vy_af = vy_ap + (lay1 + 2 * (lay2+lay3) + lay4) * 0.16666666666;

		z_tf = z_tp + (ktz1 + 2 * (ktz2+ktz3) + ktz4) * 0.16666666666;
		z_sf = z_sp + (ksz1 + 2 * (ksz2+ksz3) + ksz4) * 0.16666666666;
		z_lf = z_lp + (klz1 + 2 * (klz2+klz3) + klz4) * 0.16666666666;
		z_af = z_ap + (kaz1 + 2 * (kaz2+kaz3) + kaz4) * 0.16666666666;

		vz_tf = vz_tp + (ltz1 + 2 * (ltz2+ltz3) + ltz4) * 0.16666666666;
		vz_sf = vz_sp + (lsz1 + 2 * (lsz2+lsz3) + lsz4) * 0.16666666666;
		vz_lf = vz_lp + (llz1 + 2 * (llz2+llz3) + llz4) * 0.16666666666;
		vz_af = vz_ap + (laz1 + 2 * (laz2+laz3) + laz4) * 0.16666666666;
		fprintf(in,"%f \t %f \t %f\n",x_tf,y_tf,z_tf);
		
		x_tp =	x_tf;
		x_sp =	x_sf;
		x_lp =	x_lf;
		x_ap =	x_af;
		vx_tp =	vx_tf;
		vx_sp =	vx_sf;
		vx_lp =	vx_lf;
		vx_ap =	vx_af;
		y_tp =	y_tf;
		y_sp =	y_sf;
		y_lp =	y_lf;
		y_ap =	y_af;
		vy_tp =	vy_tf;
		vy_sp =	vy_sf;
		vy_lp =	vy_lf;
		vy_ap =	vy_af;
		z_tp =	z_tf;
		z_sp =	z_sf;
		z_lp =	z_lf;
		z_ap =	z_af;
		vz_tp =	vz_tf;
		vz_sp =	vz_sf;
		vz_lp =	vz_lf;
		vz_ap =	vz_af;
		
		
	}
		fclose(in);
}

int main(){
	rk4();
}