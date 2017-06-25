#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define HI(num) (((num) & 0x0000FF00) >> 8)
#define LO(num) ((num) & 0x000000FF)
/* run this program using the console pauser or add your own getch, system("pause") or input loop */

//https://ugurkoltuk.wordpress.com/2010/03/04/an-extreme-simple-pgm-io-api/ 

typedef struct _PGMData {
    int row;
    int col;
    int max_gray;
    int **matrix;
} PGMData;

int **allocate_dynamic_matrix(int row, int col)
{
    int **ret_val;
    int i;
 
    ret_val = (int **)malloc(sizeof(int *) * row);
    if (ret_val == NULL) {
        perror("memory allocation failure");
        exit(EXIT_FAILURE);
    }
 
    for (i = 0; i < row; ++i) {
        ret_val[i] = (int *)malloc(sizeof(int) * col);
        if (ret_val[i] == NULL) {
            perror("memory allocation failure");
            exit(EXIT_FAILURE);
        }
    }
 
    return ret_val;
}
 
void deallocate_dynamic_matrix(int **matrix, int row)
{
    int i;
 
    for (i = 0; i < row; ++i)
        free(matrix[i]);
    free(matrix);
} 
void SkipComments(FILE *fp)
{
    int ch;
    char line[100];
 
    while ((ch = fgetc(fp)) != EOF && isspace(ch))
        ;
    if (ch == '#') {
        fgets(line, sizeof(line), fp);
        SkipComments(fp);
    } else
        fseek(fp, -1, SEEK_CUR);
}
/*for reading:*/
PGMData* readPGM(const char *file_name, PGMData *data)
{
    FILE *pgmFile;
    char version[3];
    int i, j;
    int lo, hi;
 
    pgmFile = fopen(file_name, "rb");
    if (pgmFile == NULL) {
        perror("cannot open file to read");
        exit(EXIT_FAILURE);
    }
 
    fgets(version, sizeof(version), pgmFile);
    if (strcmp(version, "P5")) {
        fprintf(stderr, "Wrong file type!\n");
        exit(EXIT_FAILURE);
    }
 
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->col);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->row);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->max_gray);
    fgetc(pgmFile);
 
    data->matrix = allocate_dynamic_matrix(data->row, data->col);
    if (data->max_gray > 255)
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                hi = fgetc(pgmFile);
                lo = fgetc(pgmFile);
                data->matrix[i][j] = (hi << 8) + lo;
            }
    else
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                lo = fgetc(pgmFile);
                data->matrix[i][j] = lo;
            }
 
    fclose(pgmFile);
    return data;
 
}
 
/*and for writing*/
 
void writePGM(const char *filename, const PGMData *data)
{
    FILE *pgmFile;
    int i, j;
    int hi, lo;
 
    pgmFile = fopen(filename, "wb");
    if (pgmFile == NULL) {
        perror("cannot open file to write");
        exit(EXIT_FAILURE);
    }
 
    fprintf(pgmFile, "P5 ");
    fprintf(pgmFile, "%d %d ", data->col, data->row);
    fprintf(pgmFile, "%d ", data->max_gray);
 
    if (data->max_gray > 255) {
        for (i = 0; i < data->row; ++i) {
            for (j = 0; j < data->col; ++j) {
                hi = HI(data->matrix[i][j]);
                lo = LO(data->matrix[i][j]);
                fputc(hi, pgmFile);
                fputc(lo, pgmFile);
            }
 
        }
    } else {
        for (i = 0; i < data->row; ++i)
            for (j = 0; j < data->col; ++j) {
                lo = LO(data->matrix[i][j]);
                fputc(lo, pgmFile);
            }
    }
 
    fclose(pgmFile);
    deallocate_dynamic_matrix(data->matrix, data->row);
}

// averaging filter is applying image
void averagingFilter(const PGMData *data,int filterSize,PGMData *newData){
	int i,j,k=0,l=0,offset;
	offset=filterSize/2;
	int toplam=0;
	for(i=0;i<(data->row-filterSize);i++){
		for(j=0;j<data->col-filterSize;j++){
		 	while(k<filterSize){
				while(l<filterSize){
				    toplam=toplam+data->matrix[i+k][j+l];
						l++;		
				}
				l=0;
			k++;
		} 
 k=0;
     	 newData->matrix[i+offset][j+offset]=(toplam/(filterSize*filterSize));
		toplam=0;
	}
}
	
}

//creating gaussion kernel by filter and sigma size
void createGaussionKernel(double ** kernel,double sigma,int filterSize){
	double s=(sigma*sigma)*2.0,sum=0.0,c;
	int i,j;
 	int offset=(filterSize-1)/2;
	for(i=0;i<filterSize;i++){
	 
		for(j=0;j< filterSize;j++){
            kernel[i][j] = (exp(-((i-offset)*(i-offset)+(j-offset)*(j-offset))/s));	
		}
	}
	
	c = 1.0/ kernel[0][0];
   
 for(i=0; i< filterSize; i++)
    {
        for(j=0; j<filterSize; j++)
        {          
            kernel[i][j]*=c;
            sum +=  kernel[i][j];
        }
    }
     for(i=0; i< filterSize; i++)
    {
        for(j=0; j<filterSize; j++)
        {
            kernel[i][j]/= sum ;
        }
    }
}

//gaussing filter is applying over image
void gaussianFilter(const PGMData *data,int filterSize,PGMData *newData,double ** kernel){
	int i,j,k=0,l=0,offset;
	offset=filterSize/2;

	double toplam=0;
	for(i=0;i<(data->row-filterSize);i++){
		for(j=0;j<(data->col-filterSize);j++){
		 	while(k<filterSize){
				while(l<filterSize){
				    toplam=toplam+(data->matrix[i+k][j+l]*kernel[k][l]);
						l++;		
				}
				l=0;
			k++;
		} 
 k=0;
     	 newData->matrix[i+offset][j+offset]=(toplam);
		toplam=0;
	}
}
	
}

int main(int argc, char *argv[]) {

	int type=0,filterSize=0,i=0;
	char filename[100];
	double sigma=0.0;
	PGMData *data=(PGMData*)malloc(sizeof(PGMData));
	PGMData *newData=(PGMData*)malloc(sizeof(PGMData));
	while(1){
		printf("Select filter type: \n");
		printf("1.Averaging Filter \n");
		printf("2.Gauss Filter \n");
		printf("0. exit \n");
		scanf("%d",&type);
		
		printf("\nEnter the file name:  \n");
		scanf("%s",filename);
		printf("%s",filename);
		printf("\nEnter the filter size:  \n");
		scanf("%d",&filterSize);	
		
		i=0;
		if(type==1){

	
			readPGM(filename,data);
			newData->row=data->row;
			newData->col=data->col;
			newData->max_gray=data->max_gray;
			newData->matrix = allocate_dynamic_matrix(newData->row, newData->col); 
		    averagingFilter(data,filterSize,newData);
 		   	writePGM("filteredImage.pgm",newData);
 		   	printf("\nFiltered image created\n");
			
			
			
		}else if(type==2){

			printf("Sigma degeri girin:  \n");
			scanf("%lf",&sigma);
	        double **gaussianKernel = (double **)malloc(sizeof(double *) * filterSize);
  
		    for (i = 0; i < filterSize; ++i) {
		        gaussianKernel[i] = (double *)malloc(sizeof(double) * filterSize);
		
		    }
		 
		 	readPGM(filename,data);
			newData->row=data->row;
			newData->col=data->col;
			newData->max_gray=data->max_gray;
			newData->matrix = allocate_dynamic_matrix(newData->row, newData->col); 
		 
		  	createGaussionKernel(gaussianKernel,sigma,filterSize);
		 	gaussianFilter(data,filterSize,newData,gaussianKernel);
 	
  			writePGM("filteredImage.pgm",newData);
  		   	printf("\nFiltered image created\n");
  		   	
	   	    for (i = 0; i < filterSize; ++i) {
		    	    free(gaussianKernel[i]);
			    }
			    free(gaussianKernel);
		}else {
			
			return -1;
		}
	 
		
	}
}

