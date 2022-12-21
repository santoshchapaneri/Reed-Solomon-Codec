#include <stdio.h>
#include <stdlib.h>
#include "field.h"
#include "functions.h"

int ComputeSyndrome(GFPowRep *block, GFPowRep* Syndromes);
void ErrorLocator(GFPowRep *Syndromes, GFPowRep *ErrorLocations, GFPowRep *ErrorPoly);
void findErrorMagnitude(GFPowRep *Syndromes, GFPowRep *ErrorLocations,GFPowRep *ErrorPoly,GFPowRep *ErrMag);


void RSDecoder(char *fin, char *fout)
{
	FILE *finput, *foutput;
	int i, j, k, breakflag = 0, synerror = 0;
	unsigned char input;
	GFPowRep *Syndromes, *block, *ErrorLocations, *ErrorPoly, *ErrMag;

	// Open the files
	if((finput = fopen(fin,"rb")) == NULL)
	{
		printf("\nInput (encoded) file not found!\n");
		return;
	}
	
	if((foutput = fopen(fout,"wb")) == NULL)
	{
		printf("\nCannot create output file!\n");
		return;
	}

	//Computing Syndromes
	Syndromes = (GFPowRep*)malloc(PAR*sizeof(GFPowRep));
	block = (GFPowRep*)malloc((MSG+PAR)*sizeof(GFPowRep));
	ErrorPoly = (GFPowRep*)malloc((PAR + 1)*sizeof(GFPowRep));
	ErrorLocations = (GFPowRep*)malloc((PAR/2)*sizeof(GFPowRep));
	ErrMag = (GFPowRep*)malloc((PAR/2)*sizeof(GFPowRep));
	
//	while(1)
	{
/*		for(i=0;i<MSG+PAR;i++)
		{
			fscanf(finput,"%c",&input);
			if(feof(finput))
			{
				for(j=i-PAR,k=0;j<i;j++,k++)
				{
					block[MSG+k] = block[j];
				}
				for(j=i-PAR;j<MSG;j++)
					block[j] = (GFPowRep)46;	// '.' character
				breakflag = 1;
				break;
			}
			else
			{
				input = input%FIELD;
				block[i] = (GFPowRep)(input);
			}
		}
*/
		block[0] = 0; block[1] = 15; block[2] = 15; block[3] = 0;
		block[4] = 0; block[5] = 0; block[6] = 15; block[7] = 15;
		block[8] = 15; block[9] = 15; block[10] = 0; block[11] = 15;
		block[12] = 0; block[13] = 15; block[14] = 15; 
		//Compute Syndromes
		synerror = ComputeSyndrome(block, Syndromes);
		
		//Display Syndromes
		for(j=0;j<PAR;j++)
			printf("Syn[%d] = %d ",j,Syndromes[j]);
		printf("\n");

		if(synerror) //Error in Received Block
		{
			//Berlekamp Massey Algorithm
			ErrorLocator(Syndromes, ErrorLocations, ErrorPoly);
			
			//Forney Algorithm - determine error magnitudes
			findErrorMagnitude(Syndromes, ErrorLocations, ErrorPoly, ErrMag);
			
			//Correct the Erroneous Block
			for(j=0;j<PAR/2;j++)
			{
				if(ErrorLocations[j]!=FIELD-1)
				{
					block[ErrorLocations[j]] = GFAddition(block[ErrorLocations[j]],ErrMag[j]);
				}
			}

		}

		//Write the (corrected) block to output file
		//for(j=0;j<i-PAR;j++)
		for(j=0;j<9;j++)
		{
			fprintf(foutput,"%c",block[j]);
		}

//		if(breakflag) //End of file
//			break;
	}

	// Close the files
	fclose(finput);
	fclose(foutput);
}

int ComputeSyndrome(GFPowRep *block, GFPowRep* Syndromes)
{
	int i, j;

	//Initialize to "ZERO"
	for(j=0;j<PAR;j++)
	{
		Syndromes[j] = FIELD-1;
	}

	for(j=0;j<PAR;j++)
		Syndromes[j] = GFPolyEvaluate(block,MSG+PAR-1,(GFPowRep)(j+1));

	for(i=0;i<PAR;i++)
	{
		if(Syndromes[i]!=FIELD-1)
			return 1;//Error in recd block
	}
	return 0; //No Error
}

void ErrorLocator(GFPowRep *Syndromes, GFPowRep *ErrorLocations, GFPowRep *ErrorPoly)
{
	int i, j, k, L; 
	GFPowRep sum=FIELD-1;
	GFPowRep Delta, temp, *Tx, *PreviousEL;

	Tx = (GFPowRep*)malloc((PAR + 1)*sizeof(GFPowRep));
	PreviousEL = (GFPowRep*)malloc((PAR + 1)*sizeof(GFPowRep));

	//Initialization
	k = 0; L = 0;
	ErrorPoly[0] = 0;//alpha^0=1
	PreviousEL[0] = 0;
	Tx[0] = FIELD-1;
	for(i=1;i<(PAR)+1;i++)
	{
		ErrorPoly[i] = FIELD-1;//"ZERO"
		PreviousEL[i] = FIELD-1;
		Tx[i] = FIELD-1;
	}
	Tx[1] = 0;//alpha^0=1
	
	printf("k\tSk\tEL(0,1,2,3,4,5,6)\tDelta\tL\tTx(0,1,2,3,4,5,6)\n\t\t");

	for(k=1;k<=PAR;k++)
	{
		sum = FIELD-1;
		for(i=1;i<=L;i++)
		{
			temp = GFMultiplication(PreviousEL[i],Syndromes[k-i-1]);
			sum = GFAddition(temp,sum); 
		}

		Delta = GFAddition(Syndromes[k-1],sum);
		if(Delta != FIELD-1)
		{
			//modify connection poly
			for(j=0;j<(PAR)+1;j++)
			{
				temp = GFMultiplication(Delta,Tx[j]);
				ErrorPoly[j] = GFAddition(PreviousEL[j],temp);
			}
			
			if((2*L)<k)
			{
				L = k - L;
				temp = GFExponentiation(Delta,-1);
				for(j=0;j<(PAR/2)+1;j++)
				{
					Tx[j] = GFMultiplication(PreviousEL[j],temp);
				}
			}
		}

		for(j=(PAR);j>0;j--)
			Tx[j] = Tx[j-1];
		
		Tx[0] = FIELD-1;

		for(j=0;j<(PAR)+1;j++)
			PreviousEL[j] = ErrorPoly[j];
		
		printf("\n%d\t%d\t%d,%d,%d,%d,%d,%d,%d\t%d\t%d\t%d,%d,%d,%d,%d,%d,%d",k,Syndromes[k-1],
			ErrorPoly[0],ErrorPoly[1],ErrorPoly[2],ErrorPoly[3],ErrorPoly[4],
			ErrorPoly[5],ErrorPoly[6],Delta,L,
			Tx[0],Tx[1],Tx[2],Tx[3],Tx[4],Tx[5],Tx[6]);
	}
	printf("\n");

	for(k=0;k<PAR/2;k++)
		ErrorLocations[k] = FIELD-1;

	k = 0;
	sum = FIELD-1;

	//Find roots of ErrorPoly
	for(j=0;j<MSG+PAR;j++)
	{
		sum = GFPolyEvaluate(ErrorPoly,PAR,(GFPowRep)j);

		if(sum == FIELD-1)
		{
			ErrorLocations[k] = (GFPowRep)(j);
			k++;
		}
		else
			sum = FIELD-1;
	}
	
	for(i=0;i<k;i++) //max k = PAR
	{
		//take reciprocal
		ErrorLocations[i] = GFExponentiation(ErrorLocations[i],-1);
		printf("Error[%d]=%d\n",i,ErrorLocations[i]);
	}
}

void findErrorMagnitude(GFPowRep *Syndromes, GFPowRep *ErrorLocations,GFPowRep *ErrorPoly,GFPowRep *ErrMag)
{
	int i, j;

	GFPowRep *ErrMagPoly, *temp, *DerErrPoly;
	GFPowRep den, Xkinv;

	ErrMagPoly = (GFPowRep*)malloc((PAR+1)*sizeof(GFPowRep));
	temp = (GFPowRep*)malloc((PAR+1)*sizeof(GFPowRep));
	DerErrPoly = (GFPowRep*)malloc((PAR)*sizeof(GFPowRep));

	for(i=0;i<PAR+1;i++)
		ErrMagPoly[i] = FIELD-1;

	temp[0] = 0;
	for(i=1;i<=PAR;i++)
		temp[i] = Syndromes[i-1];
	
	for(j=0;j<PAR+1;j++)
	{
		for(i=0;i<PAR+1-j;i++)
			ErrMagPoly[j+i] = GFAddition(GFMultiplication(ErrorPoly[j],temp[i]),ErrMagPoly[j+i]);			
	}

	printf("\nLocation\tMagnitude\n");

	//Finding Derivative of ErrorPoly
	for(i=0;i<PAR;i++)
	{
		if((i%2) == 0)
			DerErrPoly[i] = ErrorPoly[i+1];
		else
			DerErrPoly[i] = FIELD-1;
	}

	for(i=0;i<PAR/2;i++)
	{
		if(ErrorLocations[i]!=FIELD-1)
		{
			Xkinv = GFExponentiation(ErrorLocations[i],-1);
			den = GFPolyEvaluate(DerErrPoly,PAR-1,Xkinv);
			den = GFExponentiation(den,-1);
			ErrMag[i] = GFPolyEvaluate(ErrMagPoly,PAR,Xkinv);
			ErrMag[i] = GFMultiplication(ErrMag[i],ErrorLocations[i]);
			ErrMag[i] = GFMultiplication(ErrMag[i],den);
		}
		else
			ErrMag[i] = FIELD-1;

		printf("   %d\t\t   %d\n",ErrorLocations[i],ErrMag[i]);
	}
}

GFPowRep GFPolyEvaluate(GFPowRep *Poly, int degree, GFPowRep arg)
{
	int i;
	GFPowRep temp, sum;

	sum = FIELD-1;

	for(i=0;i<=degree;i++)
	{
		temp = GFExponentiation(arg,i);
		temp = GFMultiplication(Poly[i],temp);
		sum = GFAddition(sum,temp);
	}

	return sum;
}