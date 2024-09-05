#include "stdafx.h"
#include<stdio.h>
#include<string.h>
#include<limits.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<float.h>
#include<Windows.h>



/*
=============================================================SPEECH REPRESENTATTION MODULE=====================================================================================================
*/

//Normalize the data
void normalize_data(char file[100]){
	//open inputfile
	FILE* fp=fopen(file,"r");
	if(fp==NULL){
		printf("Error in Opening File!\n");
		return;
	}
	int amp=0,avg=0;
	int i=0;
	int n=0;
	int min_amp=INT_MAX;
	int max_amp=INT_MIN;
	//calculate average, minimum & maximum amplitude
	while(!feof(fp)){
		fscanf(fp,"%d",&amp);
		avg+=amp;
		min_amp=(amp<min_amp)?amp:min_amp;
		max_amp=(amp>max_amp)?amp:max_amp;
		n++;
	}
	avg/=n;
	T=(n-FS)/80 + 1;
	if(T>T_) T=T_;
	//update minimum & maximum amplitude after DC Shift
	min_amp-=avg;
	max_amp-=avg;
	fseek(fp,0,SEEK_SET);
	while(!feof(fp)){
		fscanf(fp,"%d",&amp);
		if(min_amp==max_amp){
			amp=0;
		}
		else{
			//handle DC Shift
			amp-=avg;
			//normalize the data
			amp=(amp*CLIP)/((max_amp>min_amp)?max_amp:(-1)*min_amp);
			//store normalized data
			samples[i++]=amp;
		}
	}
	fclose(fp);
}

//calculate energy of frame
void calculate_energy_of_frame(int frame_no){
	int sample_start_index=frame_no*80;
	energy[frame_no]=0;
	for(int i=0;i<FS;i++){
		energy[frame_no]+=samples[i+sample_start_index]*samples[i+sample_start_index];
		energy[frame_no]/=FS;
	}
}

//Calculate Max Energy of file
long double calculate_max_energy(){
	int nf=T;
	long double max_energy=DBL_MIN;
	for(int f=0;f<nf;f++){
		if(energy[f]>max_energy){
			max_energy=energy[f];
		}
	}
	return max_energy;
}

//calculate average energy of file
long double calculate_avg_energy(){
	int nf=T;
	long double avg_energy=0.0;
	for(int f=0;f<nf;f++){
		avg_energy+=energy[f];
	}
	return avg_energy/nf;
}

//mark starting and ending of speech activity
void mark_checkpoints(){
	int nf=T;
	//Calculate energy of each frame
	for(int f=0;f<nf;f++){
		calculate_energy_of_frame(f);
	}
	//Make 10% of average energy as threshold
	long double threshold_energy=calculate_avg_energy()/10;
	//long double threshold_energy=calculate_max_energy()/10;
	int isAboveThresholdStart=1;
	int isAboveThresholdEnd=1;
	start_frame=0;
	end_frame=nf-1;
	//Find start frame where speech activity starts
	for(int f=0;f<nf-5;f++){
		for(int i=0;i<5;i++){
			isAboveThresholdStart*=(energy[f+i]>threshold_energy);
		}
		if(isAboveThresholdStart){
			start_frame=((f-5) >0)?(f-5):(0);
			break;
		}
		isAboveThresholdStart=1;
	}
	//Find end frame where speech activity ends
	for(int f=nf-1;f>4;f--){
		for(int i=0;i<5;i++){
			isAboveThresholdEnd*=(energy[f-i]>threshold_energy);
		}
		if(isAboveThresholdEnd){
			end_frame=((f+5) < nf)?(f+5):(nf-1);
			break;
		}
		isAboveThresholdEnd=1;
	}
}

//load codebook
void load_codebook(){
	FILE* fp;
	//fp=fopen("digit_codebook.csv","r");
	fp=fopen("234101053_codebook.csv","r");
	if(fp==NULL){
		printf("Error in Opening File!\n");
		return;
	}
	for(int i=1;i<=M;i++){
		for(int j=1;j<=Q;j++){
			fscanf(fp,"%Lf,",&reference[i][j]);
		}
	}
	fclose(fp);
}

//Calculate ai's using Durbin's Algo
void durbinAlgo(){
	//step-0:initialize energy
	long double E=R[0];
	long double alpha[13][13];
	for(int i=1;i<=P;i++){
		double k;
		long double numerator=R[i];
		long double alphaR=0.0;
		for(int j=1;j<=(i-1);j++){
			alphaR+=alpha[j][i-1]*R[i-j];
		}
		numerator-=alphaR;
		//step-1: calculate k
		k=numerator/E;
		//step-2: calculate alpha[i][i]
		alpha[i][i]=k;
		//step-3: calculate alpha[j][i]
		for(int j=1;j<=(i-1);j++){
			alpha[j][i]=alpha[j][i-1]-(k*alpha[i-j][i-1]);
			if(i==P){
				a[j]=alpha[j][i];
			}
		}
		//step-4: update energy
		E=(1-k*k)*E;
		if(i==P){
			a[i]=alpha[i][i];
		}
	}
}

//Calculate minimun LPC Coefficients using AutoCorrelation
void autoCorrelation(int frame_no){
	long double s[FS];
	int sample_start_index=frame_no*80;
	
	//Hamming Window Function
	for(int i=0;i<FS;i++){
		long double wn=0.54-0.46*cos((2*(22.0/7.0)*i)/(FS-1));
		s[i]=wn*samples[i+sample_start_index];
	}
	
	//Calculate R0 to R12
	for(int i=0;i<=P;i++){
		long double sum=0.0;
		for(int y=0;y<=FS-1-i;y++){
			sum+=((s[y])*(s[y+i]));
		}
		R[i]=sum;
	}

	//Apply Durbin's Algorithm to calculate ai's
	durbinAlgo();
}


//Apply Cepstral Transformation to LPC to get Cepstral Coefficient
void cepstralTransformation(){
	C[0]=2.0*(log(R[0])/log(2.0));
	for(int m=1;m<=P;m++){
		C[m]=a[m];
		for(int k=1;k<m;k++){
			C[m]+=((k*C[k]*a[m-k])/m);
		}
	}
}

//Apply raised Sine window on Cepstral Coefficients
void raisedSineWindow(){
	for(int m=1;m<=P;m++){
		//raised sine window
		long double wm=(1+(Q/2)*sin(pie*m/Q));
		C[m]*=wm;
	}
}

//Store Cepstral coefficients of each frame of file
void process_universe_file(FILE* fp, char file[]){
	//normalize data
	normalize_data(file);
	int m=0;
	int nf=T;
	//repeat procedure for frames
	for(int f=0;f<nf;f++){
		//Apply autocorrelation
		autoCorrelation(f);
		//Apply cepstral Transformation
		cepstralTransformation();
		//apply raised sine window "or" liftering
		raisedSineWindow();
		for(int i=1;i<=Q;i++){
			fprintf(fp,"%Lf,",C[i]);
		}
		fprintf(fp,"\n");
		//printf(".");
	}
}

//Generate Universe from given dataset
void generate_universe(){
	FILE* universefp;
	universefp=fopen("234101053_universe.csv","w");
	for(int d=0;d<=9;d++){
		for(int u=1;u<=TRAIN_SIZE;u++){
			char filename[40];
			_snprintf(filename,40,"234101053_dataset/234101053_E_%d_%d.txt",d,u);
			process_universe_file(universefp,filename);
		}
	}
	
}

//calculate minimium Tokhura Distance
int minTokhuraDistance(long double testC[]){
	long double minD=DBL_MAX;
	int minDi=0;
	for(int i=1;i<=M;i++){
		long double distance=0.0;
		for(int j=1;j<=Q;j++){
			distance+=(tokhuraWeight[j]*(testC[j]-reference[i][j])*(testC[j]-reference[i][j]));
		}
		if(distance<minD){
			minD=distance;
			minDi=i;
		}
	}
	return minDi;
}

//Generate Observation Sequence
void generate_observation_sequence(char file[]){
	FILE* fp=fopen("o.txt","w");
	//normalize data
	normalize_data(file);
	int m=0;
	//mark starting and ending index
	mark_checkpoints();
	T=(end_frame-start_frame+1);
	int nf=T;
	//long double avg_energy=calculate_avg_energy();
	//repeat procedure for each frames
	for(int f=start_frame;f<=end_frame;f++){
		//Apply autocorrelation
		autoCorrelation(f);
		//Apply cepstral Transformation
		cepstralTransformation();
		//apply raised sine window "or" liftering
		raisedSineWindow();
		fprintf(fp,"%d ",minTokhuraDistance(C));
	}
	fprintf(fp,"\n");
	fclose(fp);
}


/*
================================================================LBG MODULE========================================================================================================
*/
/*
	Load Universe Vector from file
	@param : file -> store filename
*/
void load_universe(char file[100]){
	//open inputfile
	FILE* fp=fopen(file,"r");
	if(fp==NULL){
		printf("Error in Opening File!\n");
		return;
	}
	
	int i=0;
	long double c;
	while(!feof(fp)){
		fscanf(fp,"%Lf,",&c);
		X[LBG_M][i]=c;
		//Ceptral coeffecient index
		i=(i+1)%12;
		//Compute Universe vector size
		if(i==0) LBG_M++;
	}
	fclose(fp);
}

/*
	Store Codebook to file
	@param : k->size of codebook
*/
void store_codebook(char file[100],int k){
	FILE* fp=fopen(file,"w");
	if(fp==NULL){
		printf("Error opening file!\n");
		return;
	}
	for(int i=0;i<k;i++){
		for(int j=0;j<12;j++){
			fprintf(fp,"%Lf,",codebook[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


/*
	Print Codebook on terminal
	@param : k->size of codebook
*/
void print_codebook(int k){
	printf("Codebook of size %d:\n",k);
	for(int i=0;i<k;i++){
		for(int j=0;j<12;j++){
			printf("%Lf\t",codebook[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


/*
	Initialize codebook with centroid of the Universe
*/
void initialize_with_centroid(){
	long double centroid[12]={0.0};
	for(int i=0;i<LBG_M;i++){
		for(int j=0;j<12;j++){
			centroid[j]+=X[i][j];
		}
	}
	for(int i=0;i<12;i++){
		centroid[i]/=LBG_M;
		codebook[0][i]=centroid[i];
	}
	print_codebook(1);
}


/*
	Calculate distance between input and codevector
	@param : x[]->input vector, y[]->codevector
	@return : long double distance between input & codevector
*/
long double calculate_distance(long double x[12], long double y[12]){
	long double distance=0.0;
	for(int i=0;i<12;i++){
		distance+=(tokhuraWeight[i+1]*(x[i]-y[i])*(x[i]-y[i]));
	}
	return distance;
}


/*
	Classification of Universe into k clusters
	@param : k->size of codebook
*/
void nearest_neighbour(int k){
	for(int i=0;i<LBG_M;i++){
		//store minimum distance between input and codebook
		long double nn=DBL_MAX;
		//store index of codevector with which input has minimum distance
		int cluster_index;
		for(int j=0;j<k;j++){
			//compute distance between input and codevector
			long double dxy=calculate_distance(X[i],codebook[j]);
			if(dxy<=nn){
				cluster_index=j;
				nn=dxy;
			}
		}
		//classification of ith input to cluster_index cluster
		cluster[i]=cluster_index;
	}
}


/*
	codevector updation
	@param : k->size of codebook
*/
void codevector_update(int k){
	long double centroid[K][12]={0.0};
	//Store number of vectors in each cluster
	int n[K]={0};
	for(int i=0;i<LBG_M;i++){
		for(int j=0;j<12;j++){
			centroid[cluster[i]][j]+=X[i][j];
		}
		n[cluster[i]]++;
	}
	//Codevector Updation as Centroid of each cluster
	for(int i=0;i<k;i++){
		for(int j=0;j<12;j++){
			codebook[i][j]=centroid[i][j]/n[i];
		}
	}
}


/*
	Calculate overall average Distortion
	@return : Distortion/ D
	D=(1/M)*sigma_1toM(d(x(n),y(n)))
*/
long double calculate_distortion(){
	long double distortion=0.0;
	for(int i=0;i<LBG_M;i++){
		distortion+=calculate_distance(X[i],codebook[cluster[i]]);
	}
	distortion/=LBG_M;
	return distortion;
}


/*
	K-Means Algorithm
	@param : k->size of codebook
*/
void KMeans(int k){
	FILE* fp=fopen("distortion.txt","a");
	if(fp==NULL){
		printf("Error pening file!\n");
		return;
	}
	//iterative index
	int m=0;
	//store previous and current D
	long double prev_D=DBL_MAX, cur_D=DBL_MAX;
	//repeat until convergence
	do{
		//Classification
		nearest_neighbour(k);
		//Iterative index update
		m++;
		//Codevector Updation
		codevector_update(k);
		prev_D=cur_D;
		//Calculate overall avg Distortion / D
		cur_D=calculate_distortion();
		printf("m=%d\t:\t",m);
		printf("Distortion:%Lf\n",cur_D);
		fprintf(fp,"%Lf\n",cur_D);
	}while((prev_D-cur_D)>DELTA);//repeat until distortion difference is >delta
	//Print Updated Codebook
	printf("Updated ");
	print_codebook(k);
	fclose(fp);
}


/*
	LBG Algorithm
*/
void LBG(){
	printf("\nLBG Algorithm:\n");
	//Start from single codebook
	int k=1;
	//Compute codevector as centroid of universe
	initialize_with_centroid();
	//repeat until desired size codebook is reached
	while(k!=K){
		//Split each codebook entry Yi to Yi(1+epsilon) & Yi(1-epsilon)
		for(int i=0;i<k;i++){
			for(int j=0;j<12;j++){
				long double Yi=codebook[i][j];
				//Yi(1+epsilon)
				codebook[i][j]=Yi-EPSILON;
				//Yi(1-epsilon)
				codebook[i+k][j]=Yi+EPSILON;
			}
		}
		//Double size of codebook
		k=k*2;
		//Call K-means with split codebook
		KMeans(k);
	}
}

//Generate Codebook From given Universe
void generate_codebook(){
	load_universe("234101053_universe.csv");
	LBG();
	store_codebook("234101053_codebook.csv",K);
}
