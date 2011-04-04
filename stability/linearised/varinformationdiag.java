//
//  build the tree from the right data!
//  
//
//  Created by Renaud Lambiotte on Tue Jul 22 2003.
//  Copyright (c) 2003 ULB. All rights reserved.
//

import java.io.*;
import java.lang.*;
import java.util.Vector;
import java.util.*;
import java.text.*;

//toGet the listing name
public class varinformationdiag
{
		
    public static void main(String args[])
    {
   				
        String record, sub="";
        int length;
		int maximum=10;
        int separator[]=new int[maximum+1];
		char SEP1=' ';
		char character[];
		int pol;
		
		DecimalFormat df = new DecimalFormat("#0.00");
		
		int[][][] link=new int[1][1][1];
		int[][] nCommunities=new int[1][1];
		
		int n0=0;
		int nNodes=-1, nPartitions=-1;
		int nTrials=50;
		
		double factor=20.0/19.0;
		int symmetry=50;
								
		//first: reads the data of the partition: number of nodes; and initialize the data structure
		try
        {
			nNodes=Integer.parseInt(args[0]);
			nPartitions=Integer.parseInt(args[1]);
									
			link=new int[nTrials][nPartitions][nNodes];
			nCommunities=new int[nTrials][nPartitions];
			
			
			for(int x=0;x<nTrials;x++)
			{
			for(int i=0;i<nNodes;i++)
			{
				for(int j=0;j<nPartitions;j++)
				{
					link[x][j][i]=-1;
					nCommunities[x][j]=0;
				}
			}
			}
			
			
			for(int file=n0;file<nPartitions;file++)
			{
				//System.out.println(file);
				int x=0;
				
				StringBuffer buffy=new StringBuffer("");
				
				buffy.append(Math.pow(factor,file-symmetry));
				buffy.append(".all.dat");
																
				FileReader read = new FileReader(buffy.toString());
				BufferedReader rw = new BufferedReader(read);  
			
				record=rw.readLine();	
				//We read through the data
				while(record!=null )
				{	
					length=record.length();
					if(length>0)
					{
						character=record.toCharArray();
               
						pol=1;  
						boolean already=false;
				
						for(int i=0; i<=length-1; i++)
						{
				
							if((character[i]==SEP1)&&pol<=maximum)
							{
								separator[pol]=i;
								pol++;
							}
			
						}
						separator[pol]=length;
				
						sub=record.substring(0, separator[1]);
						int o=Integer.parseInt(sub);
				
						sub=record.substring(separator[1]+1, separator[2]);
						int e=Integer.parseInt(sub);
					
						if(e>nCommunities[x][file])
							nCommunities[x][file]=e;
					
						link[x][file][o]=e;
					}
					else
					{
						nCommunities[x][file]++;
						x++;
					}
					
					record=rw.readLine();
				}
				
			}
		}
        catch(IOException ioe) 
        {
            System.out.println( "IO error:" + ioe );
        }
		
		//for(int oo=0;oo<nPartitions;oo++)
		//	System.out.println(oo+" "+nCommunities[oo]);
		
		for(int oo=n0;oo<nPartitions;oo++)
		{
			//System.out.println(oo);
			int pp=oo;
			double infoM=0;
			double normaM=0;
			
			for(int x1=0; x1<nTrials;x1++)
			{
				for(int x2=0;x2<x1;x2++)
				{
					double[][] matrix=new double[nCommunities[x1][oo]][nCommunities[x2][pp]];
					double[] mh=new double[nCommunities[x1][oo]];
					double[] mv=new double[nCommunities[x2][pp]];
					double norma=0;
			
					for(int i=0;i<nCommunities[x1][oo];i++)
					{
						mh[i]=0;
						for(int j=0;j<nCommunities[x2][pp];j++)
						{
							matrix[i][j]=0;
						}
					}
						
					for(int j=0;j<nCommunities[x2][pp];j++)
						mv[j]=0;
		
					for(int i=0;i<nNodes;i++)
					{
						//System.out.println(link[oo][i]+" "+nCommunities[oo]);
						//System.out.println(link[pp][i]+" "+nCommunities[pp]);
						matrix[link[x1][oo][i]][link[x2][pp][i]]++;
					}
			
					for(int i=0;i<nCommunities[x1][oo];i++)
					{
						for(int j=0;j<nCommunities[x2][pp];j++)
						{
							mh[i]+=matrix[i][j];
							mv[j]+=matrix[i][j];
							norma+=matrix[i][j];
						}
					}
					//System.out.println("norma"+norma);
			
					//for(int i=0;i<nCommunities[oo];i++)
					//	System.out.println(mh[i]);
			
					double info=0.0;
			
					for(int i=0;i<nCommunities[x1][oo];i++)
					{
						for(int j=0;j<nCommunities[x2][pp];j++)
						{
							if(matrix[i][j]>0)
								info-=2.0 * matrix[i][j] * Math.log(matrix[i][j] * norma /(mh[i] * mv[j]));
						}
					}
			
			//System.out.println(info);
			
					double infoh=0;
					double infov=0.0;
			
					for(int i=0;i<nCommunities[x1][oo];i++)
						infoh+=mh[i] * Math.log(mh[i]/norma);
				
					for(int j=0;j<nCommunities[x2][pp];j++)
						infov+=mv[j] * Math.log(mv[j]/norma);
			
					info-=(infoh + infov);
					info/=norma;
					info/=Math.log((double)(nNodes));
			
					if(info<0)
						info=0;
					
					infoM+=info;
					normaM+=1;
					
				}
			}
			infoM/=normaM;
			
			System.out.println(Math.pow(factor,oo-symmetry)+" "+infoM+" "+normaM);
		 
			//System.out.println();
		}
		
				
	}
}


