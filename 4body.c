#include <stdlib.h>
#include <math.h>
#include <stdio.h>

double accel(double x1, double x2, double x3, double x4, double G, double m1, double m2, double m3);
//double accel0(double x1, double x2, double x3, double x4, double G, double m1, double m2, double m3);

int main(int argc, char **argv){
	//rk4();
	FILE *out;
	FILE *in;

	// Indice
	int i = 0;
	int j;
	int n = 1000; 

	
	// Nombre del Archivo
	char filename1[30] = "orbitas.txt";
	char *filename2 = argv[1];

 
	// Punteros posicion 
	double x_tp, vx_tp;
	double y_tp, vy_tp;
	double z_tp, vz_tp;
	double x_ap, vx_ap;
	double y_ap, vy_ap;
	double z_ap, vz_ap;
	double x_sp, vx_sp;
	double y_sp, vy_sp;
	double z_sp, vz_sp;
	double x_lp, vx_lp;
	double y_lp, vy_lp;
	double z_lp, vz_lp;
	double x_tf, vx_tf;
	double y_tf, vy_tf;
	double z_tf, vz_tf;
	double x_af, vx_af;
	double y_af, vy_af;
	double z_af, vz_af;
	double x_sf, vx_sf;
	double y_sf, vy_sf;
	double z_sf, vz_sf;
	double x_lf, vx_lf;
	double y_lf, vy_lf;
	double z_lf, vz_lf;

	// aceleraciones y velocidades para calcular vs y xs
	double ax;
	double ay;
	double az;
	double vx;
	double vy;
	double vz;	

	// Constantes Runge-kutta (si, son tantas, para solo llamar la funcion una vez)
	double ktx1,ktx2,ktx3,ktx4,ltx1,ltx2,ltx3,ltx4;
	double ksx1,ksx2,ksx3,ksx4,lsx1,lsx2,lsx3,lsx4;
	double kax1,kax2,kax3,kax4,lax1,lax2,lax3,lax4;
	double klx1,klx2,klx3,klx4,llx1,llx2,llx3,llx4;
	double kty1,kty2,kty3,kty4,lty1,lty2,lty3,lty4;
	double ksy1,ksy2,ksy3,ksy4,lsy1,lsy2,lsy3,lsy4;
	double kay1,kay2,kay3,kay4,lay1,lay2,lay3,lay4;
	double kly1,kly2,kly3,kly4,lly1,lly2,lly3,lly4; 
	double ktz1,ktz2,ktz3,ktz4,ltz1,ltz2,ltz3,ltz4; 
	double ksz1,ksz2,ksz3,ksz4,lsz1,lsz2,lsz3,lsz4;
	double kaz1,kaz2,kaz3,kaz4,laz1,laz2,laz3,laz4;
	double klz1,klz2,klz3,klz4,llz1,llz2,llz3,llz4;

	// Step 
	double h = atof(argv[2]); //0.00149597871;
	double N = atof(argv[3]); //1000.0;

	// Constante gravitacional
	const double G = 40;

	// Masas para los cuatro cuerpos
	const double m_s = 1.0;
	const double m_t = 1E-6;
	const double m_l = 1 * pow(3.695,-8);
	const double m_a = 5 * pow(10,-11);

	// C.I. en x y v para los 4 cuerpos en las 3 dimensiones del espacio
	in = fopen(filename2, "r");
	for(i=0;i<28;i++)
    {
	    fscanf(in, "%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n", &x_tp,&x_sp,&x_ap,&x_lp,&y_tp,&y_ap,&y_sp,&y_lp,&z_tp,&z_ap,&z_sp,&z_lp,&vx_tp,&vx_sp,&vx_ap,&vx_lp,&vy_tp,&vy_ap,&vy_sp,&vy_lp,&vz_tp,&vz_ap,&vz_sp,&vz_lp);
	}
	fclose(in);

	x_tp = x_tp * 149597871;
    x_sp = x_sp * 149597871;
    x_ap = x_ap * 149597871;
    x_lp = x_lp * 149597871; 

    y_tp = y_tp * 149597871;
    y_ap = y_ap * 149597871;
    y_sp = y_sp * 149597871;
    y_lp = y_lp * 149597871;     	

    z_tp = z_tp * 149597871;
    z_ap = z_ap * 149597871;
    z_sp = z_sp * 149597871;
    z_lp = z_lp * 149597871;	

	out = fopen(filename1, "w");
	fprintf(out,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf\n",x_tp,y_tp,z_tp,x_sp,y_sp,z_sp,x_lp,y_lp,z_lp,x_ap,y_ap,z_ap);
	fclose(out);

	out = fopen(filename1, "a");
	while(i < 2*N/h)
	{
		//Runge-Kutta para x

		// Aqui calculamos k1 para xp y l1 para vxp en los cuatro cuerpos
		ax = accel(x_tp,x_sp,x_lp,x_ap,G,m_s,m_l,m_a);
		vx = vx_tp;
		ktx1 = 0.5 * h * vx;
		ltx1 = 0.5 * h * ax;
		//printf("%lf \t %lf\n",ax,vx);
		ax = accel(x_sp,x_tp,x_lp,x_ap,G,m_t,m_l,m_a);
		vx = vx_sp;
		ksx1 = 0.5 * h * vx;
		lsx1 = 0.5 * h * ax;
		ax = accel(x_ap,x_tp,x_lp,x_sp,G,m_t,m_l,m_s);
		vx = vx_ap;
		kax1 = 0.5 * h * vx;
		lax1 = 0.5 * h * ax;
		ax = accel(x_lp,x_tp,x_sp,x_ap,G,m_t,m_s,m_a);
		vx = vx_lp;
		klx1 = 0.5 * h * vx;
		llx1 = 0.5 * h * ax;

		// Aqui calculamos k2 para xp y l2 para vxp en los cuatro cuerpos
		ax = accel(x_tp+ktx1,x_sp+ksx1,x_lp+klx1,x_ap+kax1,G,m_s,m_l,m_a);
		vx = vx_tp + ltx1;
		ktx2 = 0.5 * h * vx;
		ltx2 = 0.5 * h * ax;
		ax = accel(x_sp+ksx1,x_tp+ktx1,x_lp+klx1,x_ap+kax1,G,m_t,m_l,m_a);
		vx = vx_sp + lsx1;
		ksx2 = 0.5 * h * vx;
		lsx2 = 0.5 * h * ax;
		ax = accel(x_ap+kax1,x_tp+ktx1,x_lp+klx1,x_sp+ksx1,G,m_t,m_l,m_s);
		vx = vx_ap + lax1;
		kax2 = 0.5 * h * vx;
		lax2 = 0.5 * h * ax;
		ax = accel(x_lp+klx1,x_tp+ktx1,x_sp+ksx1,x_ap+kax1,G,m_t,m_s,m_a);
		vx = vx_lp + llx1;
		klx2 = 0.5 * h * vx;
		llx2 = 0.5 * h * ax;

		// Aqui calculamos k3 para xp y l3 para vxp en los cuatro cuerpos
		ax = accel(x_tp+ktx2,x_sp+ksx2,x_lp+klx2,x_ap+kax2,G,m_s,m_l,m_a);
		vx = vx_tp + ltx2;
		ktx3 = h * vx;
		ltx3 = h * ax;
		ax = accel(x_sp+ksx2,x_tp+ktx2,x_lp+klx2,x_ap+kax2,G,m_t,m_l,m_a);
		vx = vx_sp + lsx2;
		ksx3 = h * vx;
		lsx3 = h * ax;
		ax = accel(x_ap+kax2,x_tp+ktx2,x_lp+klx2,x_sp+ksx2,G,m_t,m_l,m_s);
		vx = vx_ap + lax2;
		kax3 = h * vx;
		lax3 = h * ax;
		ax = accel(x_lp+klx2,x_tp+ktx2,x_sp+ksx2,x_ap+kax2,G,m_t,m_s,m_a);
		vx = vx_lp + llx2;
		klx3 = h * vx;
		llx3 = h * ax;

		// Aqui calculamos k4 para xp y l4 para vxp en los cuatro cuerpos
		ax = accel(x_tp+ktx3,x_sp+ksx3,x_lp+klx3,x_ap+kax3,G,m_s,m_l,m_a);
		vx = vx_tp + ltx3;
		ktx4 = 0.5 * h * vx;
		ltx4 = 0.5 * h * ax;
		ax = accel(x_sp+ksx3,x_tp+ktx3,x_lp+klx3,x_ap+kax3,G,m_t,m_l,m_a);
		vx = vx_sp + lsx3;
		ksx4 = 0.5 * h * vx;
		lsx4 = 0.5 * h * ax;
		ax = accel(x_ap+kax3,x_tp+ktx3,x_lp+klx3,x_sp+ksx3,G,m_t,m_l,m_s);
		vx = vx_ap + lax3;
		kax4 = 0.5 * h * vx;
		lax4 = 0.5 * h * ax;
		ax = accel(x_lp+klx3,x_tp+ktx3,x_sp+ksx3,x_ap+kax3,G,m_t,m_s,m_a);
		vx = vx_lp + llx3;
		klx4 = 0.5 * h * vx;
		llx4 = 0.5 * h * ax;

		//Runge-Kutta para y

		// Aqui calculamos k1 para yp y l1 para vyp en los cuatro cuerpos
		ay = accel(y_tp,y_sp,y_lp,y_ap,G,m_s,m_l,m_a);
		vy = vy_tp;
		kty1 = 0.5 * h * vy;
		lty1 = 0.5 * h * ay;
		ay = accel(y_sp,y_tp,y_lp,y_ap,G,m_t,m_l,m_a);
		vy = vy_sp;
		ksy1 = 0.5 * h * vy;
		lsy1 = 0.5 * h * ay;
		ay = accel(y_ap,y_tp,y_lp,y_sp,G,m_t,m_l,m_s);
		vy = vy_ap;
		kay1 = 0.5 * h * vy;
		lay1 = 0.5 * h * ay;
		ay = accel(y_lp,y_tp,y_sp,y_ap,G,m_t,m_s,m_a);
		vy = vy_lp;
		kly1 = 0.5 * h * vy;
		lly1 = 0.5 * h * ay;

		// Aqui calculamos k2 para yp y l2 para vyp en los cuatro cuerpos
		ay = accel(y_tp+kty1,y_sp+ksy1,y_lp+kly1,y_ap+kay1,G,m_s,m_l,m_a);
		vy = vy_tp + lty1;
		kty2 = 0.5 * h * vy;
		lty2 = 0.5 * h * ay;
		ay = accel(y_sp+ksy1,y_tp+kty1,y_lp+kly1,y_ap+kay1,G,m_t,m_l,m_a);
		vy = vy_sp + lsy1;
		ksy2 = 0.5 * h * vy;
		lsy2 = 0.5 * h * ay;
		ay = accel(y_ap+kay1,y_tp+kty1,y_lp+kly1,y_sp+ksy1,G,m_t,m_l,m_s);
		vy = vy_ap + lay1;
		kay2 = 0.5 * h * vy;
		lay2 = 0.5 * h * ay;
		ay = accel(y_lp+kly1,y_tp+kty1,y_sp+ksy1,y_ap+kay1,G,m_t,m_s,m_a);
		vy = vy_lp + lly1;
		kly2 = 0.5 * h * vy;
		lly2 = 0.5 * h * ay;

		// Aqui calculamos k3 para yp y l3 para vyp en los cuatro cuerpos
		ay = accel(y_tp+kty2,y_sp+ksy2,y_lp+kly2,y_ap+kay2,G,m_s,m_l,m_a);
		vy = vy_tp + lty2;
		kty3 = h * vy;
		lty3 = h * ay;
		ay = accel(y_sp+ksy2,y_tp+kty2,y_lp+kly2,y_ap+kay2,G,m_t,m_l,m_a);
		vy = vy_sp + lsy2;
		ksy3 = h * vy;
		lsy3 = h * ay;
		ay = accel(y_ap+kay2,y_tp+kty2,y_lp+kly2,y_sp+ksy2,G,m_t,m_l,m_s);
		vy = vy_ap + lay2;
		kay3 = h * vy;
		lay3 = h * ay;
		ay = accel(y_lp+kly2,y_tp+kty2,y_sp+ksy2,y_ap+kay2,G,m_t,m_s,m_a);
		vy = vy_lp + lly2;
		kly3 = h * vy;
		lly3 = h * ay;

		// Aqui calculamos k4 para yp y l4 para vyp en los cuatro cuerpos
		ay = accel(y_tp+kty3,y_sp+ksy3,y_lp+kly3,y_ap+kay3,G,m_s,m_l,m_a);
		vy = vy_tp + lty3;
		kty4 = 0.5 * h * vy;
		lty4 = 0.5 * h * ay;
		ay = accel(y_sp+ksy3,y_tp+kty3,y_lp+kly3,y_ap+kay3,G,m_t,m_l,m_a);
		vy = vy_sp + lsy3;
		ksy4 = 0.5 * h * vy;
		lsy4 = 0.5 * h * ay;
		ay = accel(y_ap+kay3,y_tp+kty3,y_lp+kly3,y_sp+ksy3,G,m_t,m_l,m_s);
		vy = vy_ap + lay3;
		kay4 = 0.5 * h * vy;
		lay4 = 0.5 * h * ay;
		ay = accel(y_lp+kly3,y_tp+kty3,y_sp+ksy3,y_ap+kay3,G,m_t,m_s,m_a);
		vy = vy_lp + lly3;
		kly4 = 0.5 * h * vy;
		lly4 = 0.5 * h * ay;

		//Runge Kutta para z

		// Aqui calculamos k1 para zp y l1 para vzp en los cuatro cuerpos
		az = accel(z_tp,z_sp,z_lp,z_ap,G,m_s,m_l,m_a);
		vz = vz_tp;
		ktz1 = 0.5 * h * vz;
		ltz1 = 0.5 * h * az;
		az = accel(z_sp,z_tp,z_lp,z_ap,G,m_t,m_l,m_a);
		vz = vz_sp;
		ksz1 = 0.5 * h * vz;
		lsz1 = 0.5 * h * az;
		az = accel(z_ap,z_tp,z_lp,z_sp,G,m_t,m_l,m_s);
		vz = vz_ap;
		kaz1 = 0.5 * h * vz;
		laz1 = 0.5 * h * az;
		az = accel(z_lp,z_tp,z_sp,z_ap,G,m_t,m_s,m_a);
		vz = vz_lp;
		klz1 = 0.5 * h * vz;
		llz1 = 0.5 * h * az;

		// Aqui calculamos k2 para zp y l2 para vzp en los cuatro cuerpos
		az = accel(z_tp+ktz1,z_sp+ksz1,z_lp+klz1,z_ap+kaz1,G,m_s,m_l,m_a);
		vz = vz_tp + ltz1;
		ktz2 = 0.5 * h * vz;
		ltz2 = 0.5 * h * az;
		az = accel(z_sp+ksz1,z_tp+ktz1,z_lp+klz1,z_ap+kaz1,G,m_t,m_l,m_a);
		vz = vz_sp + lsz1;
		ksz2 = 0.5 * h * vz;
		lsz2 = 0.5 * h * az;
		az = accel(z_ap+kaz1,z_tp+ktz1,z_lp+klz1,z_sp+ksz1,G,m_t,m_l,m_s);
		vz = vz_ap + laz1;
		kaz2 = 0.5 * h * vz;
		laz2 = 0.5 * h * az;
		az = accel(z_lp+klz1,z_tp+ktz1,z_sp+ksz1,z_ap+kaz1,G,m_t,m_s,m_a);
		vz = vz_lp + llz1;
		klz2 = 0.5 * h * vz;
		llz2 = 0.5 * h * az;

		// Aqui calculamos k3 para zp y l3 para vzp en los cuatro cuerpos
		az = accel(z_tp+ktz2,z_sp+ksz2,z_lp+klz2,z_ap+kaz2,G,m_s,m_l,m_a);
		vz = vz_tp + ltz2;
		ktz3 = h * vz;
		ltz3 = h * az;
		az = accel(z_sp+ksz2,z_tp+ktz2,z_lp+klz2,z_ap+kaz2,G,m_t,m_l,m_a);
		vz = vz_sp + lsz2;
		ksz3 = h * vz;
		lsz3 = h * az;
		az = accel(z_ap+kaz2,z_tp+ktz2,z_lp+klz2,z_sp+ksz2,G,m_t,m_l,m_s);
		vz = vz_ap + laz2;
		kaz3 = h * vz;
		laz3 = h * az;
		az = accel(z_lp+klz2,z_tp+ktz2,z_sp+ksz2,z_ap+kaz2,G,m_t,m_s,m_a);
		vz = vz_lp + llz2;
		klz3 = h * vz;
		llz3 = h * az;

		// Aqui calculamos k4 para zp y l4 para vzp en los cuatro cuerpos
		az = accel(z_tp+ktz3,z_sp+ksz3,z_lp+klz3,z_ap+kaz3,G,m_s,m_l,m_a);
		vz = vz_tp + ltz3;
		ktz4 = 0.5 * h * vz;
		ltz4 = 0.5 * h * az;
		az = accel(z_sp+ksz3,z_tp+ktz3,z_lp+klz3,z_ap+kaz3,G,m_t,m_l,m_a);
		vz = vz_sp + lsz3;
		ksz4 = 0.5 * h * vz;
		lsz4 = 0.5 * h * az;
		az = accel(z_ap+kaz3,z_tp+ktz3,z_lp+klz3,z_sp+ksz3,G,m_t,m_l,m_s);
		vz = vz_ap + laz3;
		kaz4 = 0.5 * h * vz;
		laz4 = 0.5 * h * az;
		az = accel(z_lp+klz3,z_tp+ktz3,z_sp+ksz3,z_ap+kaz3,G,m_t,m_s,m_a);
		vz = vz_lp + llz3;
		klz4 = 0.5 * h * vz;
		llz4 = 0.5 * h * az;

		// actualizacion de los vectores con xp hallamos xf 
		// clarificacion, si ayuda, de alguna manera para nosotros xp es present y xf es future, vivimos en el pasado. 

		x_tf = x_tp + h*(ktx1 + 2 * (ktx2+ktx3) + ktx4) * 0.16666666666;
		x_sf = x_sp + h*(ksx1 + 2 * (ksx2+ksx3) + ksx4) * 0.16666666666;
		x_lf = x_lp + h*(klx1 + 2 * (klx2+klx3) + klx4) * 0.16666666666;
		x_af = x_ap + h*(kax1 + 2 * (kax2+kax3) + kax4) * 0.16666666666;

		vx_tf = vx_tp + h*(ltx1 + 2 * (ltx2+ltx3) + ltx4) * 0.16666666666;
		vx_sf = vx_sp + h*(lsx1 + 2 * (lsx2+lsx3) + lsx4) * 0.16666666666;
		vx_lf = vx_lp + h*(llx1 + 2 * (llx2+llx3) + llx4) * 0.16666666666;
		vx_af = vx_ap + h*(lax1 + 2 * (lax2+lax3) + lax4) * 0.16666666666;

		y_tf = y_tp + h*(kty1 + 2 * (kty2+kty3) + kty4) * 0.16666666666;
		y_sf = y_sp + h*(ksy1 + 2 * (ksy2+ksy3) + ksy4) * 0.16666666666;
		y_lf = y_lp + h*(kly1 + 2 * (kly2+kly3) + kly4) * 0.16666666666;
		y_af = y_ap + h*(kay1 + 2 * (kay2+kay3) + kay4) * 0.16666666666;

		vy_tf = vy_tp + h*(lty1 + 2 * (lty2+lty3) + lty4) * 0.16666666666;
		vy_sf = vy_sp + h*(lsy1 + 2 * (lsy2+lsy3) + lsy4) * 0.16666666666;
		vy_lf = vy_lp + h*(lly1 + 2 * (lly2+lly3) + lly4) * 0.16666666666;
		vy_af = vy_ap + h*(lay1 + 2 * (lay2+lay3) + lay4) * 0.16666666666;

		z_tf = z_tp + h*(ktz1 + 2 * (ktz2+ktz3) + ktz4) * 0.16666666666;
		z_sf = z_sp + h*(ksz1 + 2 * (ksz2+ksz3) + ksz4) * 0.16666666666;
		z_lf = z_lp + h*(klz1 + 2 * (klz2+klz3) + klz4) * 0.16666666666;
		z_af = z_ap + h*(kaz1 + 2 * (kaz2+kaz3) + kaz4) * 0.16666666666;

		vz_tf = vz_tp + h*(ltz1 + 2 * (ltz2+ltz3) + ltz4) * 0.16666666666;
		vz_sf = vz_sp + h*(lsz1 + 2 * (lsz2+lsz3) + lsz4) * 0.16666666666;
		vz_lf = vz_lp + h*(llz1 + 2 * (llz2+llz3) + llz4) * 0.16666666666;
		vz_af = vz_ap + h*(laz1 + 2 * (laz2+laz3) + laz4) * 0.16666666666;

		x_tp = x_tf;
		x_sp = x_sf;
		x_lp = x_lf;
		x_ap = x_af;
		vx_tp =	vx_tf;
		vx_sp =	vx_sf;
		vx_lp =	vx_lf;
		vx_ap =	vx_af;
		y_tp = y_tf;
		y_sp = y_sf;
		y_lp = y_lf;
		y_ap = y_af;
		vy_tp =	vy_tf;
		vy_sp =	vy_sf;
		vy_lp =	vy_lf;
		vy_ap =	vy_af;
		z_tp = z_tf;
		z_sp = z_sf;
		z_lp = z_lf;
		z_ap = z_af;
		vz_tp =	vz_tf;
		vz_sp =	vz_sf;
		vz_lp =	vz_lf;
		vz_ap =	vz_af;

		fprintf(out,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf\n",x_tp,y_tp,z_tp,x_sp,y_sp,z_sp,x_lp,y_lp,z_lp,x_ap,y_ap,z_ap);

		i = i + 1 ;
		
	}
		fclose(out);
}

double accel(double x1, double x2, double x3, double x4, double G, double m1, double m2, double m3)
{
	
  	return (G * (((m1 * (x2 - x1))) + ((m2 * (x3 - x1))) + ((m3 * (x4 - x1)))));
  	

}
