/*============================================================*/
/*==================== freadHNB() ============================*/
/*============================================================*/
/* read HNBody orbital elements                               */
/*

   Richard E. Zeebe
   School of Ocean and Earth 
   Science and Technology
   University of Hawaii at Manoa
   1000 Pope Road, MSB 504
   Honolulu, HI 96822, USA
   email: zeebe@soest.hawaii.edu

   updates:

   11/20/21 appended "/" to dir: strcat(fpstr,"/") 
   12/29/20 new

   orbital elements:
   [1] t   : time
   [2] aa  : SemiMajorAxis
   [3] ee  : Eccentricity
   [4] inc : Inclination
   [5] lph : LongPerihelion
   [6] lan : LongAscendNode
   [7] arp : ArgPerihelion
   [8] mna : MeanAnomaly
  
 */
/*============================================================*/

void freadHNB(char *dir, char *foo, int *ls, double *t, double *aa,
              double *ee, double *inc, double *lph, double *lan, 
              double *arp, double *mna)
{
 int i,imax,chk=8;
 double **y;
 char *fpstr=dir,mssg[BUFSIZ],str[BUFSIZ];
 FILE *fp;

 y = dmatrix(1,NS,1,MS);
	
 /* open HNBody file */
 strcat(fpstr,"/");
 strcat(fpstr,foo);
 fp = fopen(fpstr,"r");
 if(fp == NULL){ 
    sprintf(mssg,"freadHNB(): Can't open HNBody file: '%s'",fpstr);
	ferrx(mssg);
 } 
 else{
    printf("\n@ Reading HNBody file: '%s'\n\n",fpstr);
 }

 /* use fscanf(). slightly slower than fread() at once */
 /* skip 17 header lines */
 for(i=1;i<=17;i++){	 
   fgets(str,BUFSIZ,fp);
   //printf("%s\n",str);
 }
 i = 0;
 while(chk == 8){
  i++;
  chk = fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf\n", \
    &y[1][i],&y[2][i],&y[3][i],&y[4][i],&y[5][i],&y[6][i],&y[7][i],&y[8][i]);
 }
 imax = i-1; /* -1 for EOF */

 if(chk != EOF){ 
    sprintf(mssg,"freadHNB(): Error while reading HNBody file: '%s'",fpstr);
	printf("fscanf() check value = %d \n",chk);
	ferrx(mssg);
 } 	

 /* close file */
 fclose(fp);

 /* extract variables from y */
 *ls = imax;
 for(i=1;i<=imax;i++){
     t[i] = y[1][i];
    aa[i] = y[2][i];
    ee[i] = y[3][i];
   inc[i] = y[4][i];
   lph[i] = y[5][i];
   lan[i] = y[6][i];
   arp[i] = y[7][i];
   mna[i] = y[8][i];
 }

 if(t[imax] > 0.0){ 
    sprintf(mssg,"freadHNB(): expecting tend<0 but tend = %f",t[imax]);
	ferrx(mssg);
 }
	
 free_dmatrix(y,1,NS,1,MS);	
	
}