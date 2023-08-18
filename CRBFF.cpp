#include<sys/time.h>
#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<time.h>
#include<conio.h>
#include<math.h>
#define max 10000
struct machine {
	int id,c,r,b,rank,loc;
	float pow_eff,l,h,cu,ru,bu,a,v,p_min,p_max;
};
struct machine vm[max],pm[max];
struct timeval start, end;
int i, j, x, y, nvm, npm, r, nrun=30;
void FF(struct machine [],struct machine [],int,int);
void BF(struct machine [],struct machine [],int,int);
void FFD(struct machine [],struct machine [],int,int);
void BFD(struct machine [],struct machine [],int,int);
void GRVMP(struct machine [],struct machine [],int,int);
void CRBFF(struct machine [],struct machine [],int,int);
void vm_decreasing_order(struct machine [], int);
void pm_decreasing_order_powerEfficiency(struct machine [], int);
void machines(int, int);
main()
{
	int i, j, code, nvm, npm, x, y;
	while(1) {
	//	gettint the number of pm
//		printf("\nEnter number of PMs: ");
//		scanf("%d", &npm);
		npm = 1000;
	//	gettint the number of vm
		printf("Enter number of VMs: ");
		scanf("%d", &nvm);
	//	selecting algorithms	
		while (1) {	
			printf("\n1-FF");
			printf("\n2-FFD");
			printf("\n3-BF");
			printf("\n4-BFD");
			printf("\n5-GRVMP");
			printf("\n6-CRBFF");
			printf("\n7-Start New Iteration\n"); 
			printf("\nEnter Your Choice: ");
			scanf("%d", &code);
			switch(code)
			{
				case 1: FF(pm,vm,npm,nvm);
						break;
				case 2: FFD(pm,vm,npm,nvm);
						break;
				case 3: BF(pm,vm,npm,nvm);
						break;
				case 4: BFD(pm,vm,npm,nvm);
						break;
				case 5: GRVMP(pm,vm,npm,nvm);
						break;
				case 6: CRBFF(pm,vm,npm,nvm);
						break;
				case 7: break; 	
			}
			if (code == 7)
				break;
		}
	}
}

//////////////////////////////  First Fit  //////////////////////////////
void FF(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int max_unaloc_vm=0, min_unaloc_vm=10000, sum_unaloc_vm=0; // number of unallocated vms
	int max_ac_pm=0, min_ac_pm=10000, sum_ac_pm=0;  //number of active pms
	float max_cu=0.0, min_cu=1.0, tot_cu=0.0;  //total cpu utilization
	float max_ru=0.0, min_ru=1.0, tot_ru=0.0;  //total ram utilization
	float max_rw=0.0, min_rw=1000.0, tot_rw=0.0;  //total rw
	float max_pc=0.0, min_pc=1000000.0, tot_pc=0.0;  //total pc
	float tmp_time=0;
	
	for (r=1; r<=nrun; r++) {
		int i, j, idle=0, unaloc_vm=0, ac_pm=0;
		float sum_cu=0, sum_ru=0, rw=0, pc=0, ave, var=0;
		struct machine pm_org[npm], vm_copy[nvm];
		machines(npm, nvm);
		for (i=0; i<npm; i++) 
			pm_org[i] = pm[i];
		for (i=0; i<nvm; i++) 
			vm_copy[i] = vm[i];
		
		gettimeofday(&start, NULL);
	
		for (i=0; i<nvm; i++)
			for (j=0; j<npm; j++) {
				if (vm_copy[i].c <= pm_org[j].c && vm_copy[i].r <= pm_org[j].r){
					pm_org[j].c -= vm_copy[i].c;
					pm_org[j].r -= vm_copy[i].r;
					vm_copy[i].c = 0;
					break;
				}
			}
			
		gettimeofday(&end, NULL);
	    //printf("\nElapsed time: %ld microsecons\n", (end.tv_usec - start.tv_usec));
		tmp_time += (end.tv_usec - start.tv_usec);
		
		for (i=0; i<npm; i++) {
			pm_org[i].cu = (float)(pm[i].c-pm_org[i].c)/pm[i].c;
			pm_org[i].ru = (float)(pm[i].r-pm_org[i].r)/pm[i].r;
			if (pm_org[i].c == pm[i].c)
				idle++;
			else
				rw += (fabs(pm_org[i].cu-pm_org[i].ru)+0.0001)/(pm_org[i].cu+pm_org[i].ru);
		}
		for (i=0; i<npm; i++)
			if (pm_org[i].cu > 0) {
				sum_cu += pm_org[i].cu;
				sum_ru += pm_org[i].ru;
				pc += (pm_org[i].p_min+(pm_org[i].p_max-pm_org[i].p_min)*pm_org[i].cu);
				ave = (pm_org[i].cu + pm_org[i].ru)/2;
				var += 0.5*((pm_org[i].cu - ave,2)+pow(pm_org[i].ru - ave,2));
			}
		for (i=0; i<nvm; i++) 
			if (vm_copy[i].c != 0)
				unaloc_vm++;
		ac_pm = npm - idle;

		// number of active pms	
		sum_ac_pm += ac_pm;
		if (ac_pm > max_ac_pm)
			max_ac_pm = ac_pm;
		if (ac_pm < min_ac_pm)
			min_ac_pm = ac_pm;
		// resource wastage
		tot_rw += rw;
		if (rw > max_rw)
			max_rw = rw;
		if (rw < min_rw)
			min_rw = rw;
		// power consumption 
		tot_pc += pc;
		if (pc > max_pc)
			max_pc = pc;
		if (pc < min_pc)
			min_pc = pc;
	}
	
	printf("\nResults:\n");
	printf("Algorithm: FF, Run=%d\n", nrun);
	printf("Average number of PMs used: (max=%.2f , min=%.2f , ave=%.2f)\n",max_ac_pm-(float)sum_ac_pm/nrun, (float)sum_ac_pm/nrun-min_ac_pm, (float)sum_ac_pm/nrun); 
	printf("Average Resource Wastage: (max=%.2f , min=%.2f , ave=%.2f)\n",max_rw-tot_rw/nrun, tot_rw/nrun-min_rw, tot_rw/nrun);
	printf("Average Power Consumption: (max=%.2f , min=%.2f , ave=%.2f)\n\n",max_pc-tot_pc/nrun, tot_pc/nrun-min_pc, tot_pc/nrun);
}	

//////////////////////////////  First Fit Decreasing  //////////////////////////////
void FFD(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int max_unaloc_vm=0, min_unaloc_vm=10000, sum_unaloc_vm=0; // number of unallocated vms
	int max_ac_pm=0, min_ac_pm=10000, sum_ac_pm=0;  //number of active pms
	float max_cu=0.0, min_cu=1.0, tot_cu=0.0;  //total cpu utilization
	float max_ru=0.0, min_ru=1.0, tot_ru=0.0;  //total ram utilization
	float max_rw=0.0, min_rw=1000.0, tot_rw=0.0;  //total rw
	float max_pc=0.0, min_pc=1000000.0, tot_pc=0.0;  //total pc
	float tmp_time=0;
	
	for (r=1; r<=nrun; r++) {
		int i, j, idle=0, unaloc_vm=0, ac_pm=0;
		float sum_cu=0, sum_ru=0, rw=0, pc=0, ave, var=0;
		struct machine pm_org[npm], pm_copy[npm], vm_copy[nvm];
		machines(npm, nvm);
		for (i=0; i<npm; i++) 
			pm_org[i] = pm[i];
		for (i=0; i<nvm; i++) 
			vm_copy[i] = vm[i];
			
		gettimeofday(&start, NULL);
		
		vm_decreasing_order(vm_copy, nvm);
	
		for (i=0; i<nvm; i++)
			for (j=0; j<npm; j++) {
				if (vm_copy[i].c <= pm_org[j].c && vm_copy[i].r <= pm_org[j].r){
					pm_org[j].c -= vm_copy[i].c;
					pm_org[j].r -= vm_copy[i].r;
					vm_copy[i].c = 0;
					break;
				}
			}
			
		gettimeofday(&end, NULL);
	    //printf("\nElapsed time: %ld microsecons\n", (end.tv_usec - start.tv_usec));
		tmp_time += (end.tv_usec - start.tv_usec);
		
		for (i=0; i<npm; i++) {
			pm_org[i].cu = (float)(pm[i].c-pm_org[i].c)/pm[i].c;
			pm_org[i].ru = (float)(pm[i].r-pm_org[i].r)/pm[i].r;
			if (pm_org[i].c == pm[i].c)
				idle++;
			else
				rw += (fabs(pm_org[i].cu-pm_org[i].ru)+0.0001)/(pm_org[i].cu+pm_org[i].ru);
		}
		for (i=0; i<npm; i++)
			if (pm_org[i].cu > 0) {
				sum_cu += pm_org[i].cu;
				sum_ru += pm_org[i].ru;
				pc += (pm_org[i].p_min+(pm_org[i].p_max-pm_org[i].p_min)*pm_org[i].cu);
				ave = (pm_org[i].cu + pm_org[i].ru)/2;
				var += 0.5*((pm_org[i].cu - ave,2)+pow(pm_org[i].ru - ave,2));
			}
		for (i=0; i<nvm; i++) 
			if (vm_copy[i].c != 0)
				unaloc_vm++;
		ac_pm = npm - idle;

		// number of active pms	
		sum_ac_pm += ac_pm;
		if (ac_pm > max_ac_pm)
			max_ac_pm = ac_pm;
		if (ac_pm < min_ac_pm)
			min_ac_pm = ac_pm;
		// resource wastage
		tot_rw += rw;
		if (rw > max_rw)
			max_rw = rw;
		if (rw < min_rw)
			min_rw = rw;
		// power consumption 
		tot_pc += pc;
		if (pc > max_pc)
			max_pc = pc;
		if (pc < min_pc)
			min_pc = pc;
	}
	
	printf("\nResults:\n");
	printf("Algorithm: FFD, Run=%d\n", nrun);
	printf("Average number of PMs used: (max=%.2f , min=%.2f , ave=%.2f)\n",max_ac_pm-(float)sum_ac_pm/nrun, (float)sum_ac_pm/nrun-min_ac_pm, (float)sum_ac_pm/nrun); 
	printf("Average Resource Wastage: (max=%.2f , min=%.2f , ave=%.2f)\n",max_rw-tot_rw/nrun, tot_rw/nrun-min_rw, tot_rw/nrun);
	printf("Average Power Consumption: (max=%.2f , min=%.2f , ave=%.2f)\n\n",max_pc-tot_pc/nrun, tot_pc/nrun-min_pc, tot_pc/nrun);;
}	

////////////////////////////// Best Fit //////////////////////////////
void BF(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int max_unaloc_vm=0, min_unaloc_vm=10000, sum_unaloc_vm=0; // number of unallocated vms
	int max_ac_pm=0, min_ac_pm=10000, sum_ac_pm=0;  //number of active pms
	float max_cu=0.0, min_cu=1.0, tot_cu=0.0;  //total cpu utilization
	float max_ru=0.0, min_ru=1.0, tot_ru=0.0;  //total ram utilization
	float max_rw=0.0, min_rw=1000.0, tot_rw=0.0;  //total rw
	float max_pc=0.0, min_pc=1000000.0, tot_pc=0.0;  //total pc
	float tmp_time=0;
	
	for (r=1; r<=nrun; r++) {
		int i, j, idle=0, unaloc_vm=0, ac_pm=0;
		float sum_cu=0, sum_ru=0, rw=0, pc=0, ave, var=0;
		struct machine pm_org[npm], pm_copy[npm], vm_copy[nvm];
		machines(npm, nvm);
		for (i=0; i<npm; i++) 
			pm_org[i] = pm_copy[i]= pm[i];
		for (i=0; i<nvm; i++) 
			vm_copy[i] = vm[i];
		
		gettimeofday(&start, NULL);		
		
		int min = 100000, rem, index;
		for (i=0; i<nvm; i++) {
			min = 100000;
			for (j=0; j<npm; j++) 
				if (vm_copy[i].c <= pm_org[j].c && vm_copy[i].r <= pm_org[j].r) {
					rem = pm_org[j].c - vm_copy[i].c;
					if (rem <min) {
						index = j;
						min = rem;
					}
				}
			pm_org[index].c -= vm_copy[i].c;
			pm_org[index].r -= vm_copy[i].r;
			vm_copy[i].c = 0;
		}
		
		gettimeofday(&end, NULL);
	    //printf("\nElapsed time: %ld microsecons\n", (end.tv_usec - start.tv_usec));
	    tmp_time += (end.tv_usec - start.tv_usec);
		
		for (i=0; i<npm; i++) {
			pm_org[i].cu = (float)(pm_copy[i].c-pm_org[i].c)/pm_copy[i].c;
			pm_org[i].ru = (float)(pm_copy[i].r-pm_org[i].r)/pm_copy[i].r;
			if (pm_org[i].c == pm_copy[i].c)
				idle++;
			else
				rw += (fabs(pm_org[i].cu-pm_org[i].ru)+0.0001)/(pm_org[i].cu+pm_org[i].ru);
		}
		for (i=0; i<npm; i++)
			if (pm_org[i].cu>0) {
				sum_cu += pm_org[i].cu;
				sum_ru += pm_org[i].ru;
				pc += (pm_org[i].p_min+(pm_org[i].p_max-pm_org[i].p_min)*pm_org[i].cu);
				ave = (pm_org[i].cu + pm_org[i].ru)/2;
				var += 0.5*((pm_org[i].cu - ave,2)+pow(pm_org[i].ru - ave,2));
			}
		for (i=0; i<nvm; i++) 
			if (vm_copy[i].c != 0)
				unaloc_vm++;
		ac_pm = npm - idle;
		
		// number of active pms	
		sum_ac_pm += ac_pm;
		if (ac_pm > max_ac_pm)
			max_ac_pm = ac_pm;
		if (ac_pm < min_ac_pm)
			min_ac_pm = ac_pm;
		// resource wastage
		tot_rw += rw;
		if (rw > max_rw)
			max_rw = rw;
		if (rw < min_rw)
			min_rw = rw;
		// power consumption 
		tot_pc += pc;
		if (pc > max_pc)
			max_pc = pc;
		if (pc < min_pc)
			min_pc = pc;
	}
	
	printf("\nResults:\n");
	printf("Algorithm: BF, Run=%d\n", nrun);
	printf("Average number of PMs used: (max=%.2f , min=%.2f , ave=%.2f)\n",max_ac_pm-(float)sum_ac_pm/nrun, (float)sum_ac_pm/nrun-min_ac_pm, (float)sum_ac_pm/nrun); 
	printf("Average Resource Wastage: (max=%.2f , min=%.2f , ave=%.2f)\n",max_rw-tot_rw/nrun, tot_rw/nrun-min_rw, tot_rw/nrun);
	printf("Average Power Consumption: (max=%.2f , min=%.2f , ave=%.2f)\n\n",max_pc-tot_pc/nrun, tot_pc/nrun-min_pc, tot_pc/nrun);;
}

////////////////////////////// Best Fit Decreasing //////////////////////////////
void BFD(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int max_unaloc_vm=0, min_unaloc_vm=10000, sum_unaloc_vm=0; // number of unallocated vms
	int max_ac_pm=0, min_ac_pm=10000, sum_ac_pm=0;  //number of active pms
	float max_cu=0.0, min_cu=1.0, tot_cu=0.0;  //total cpu utilization
	float max_ru=0.0, min_ru=1.0, tot_ru=0.0;  //total ram utilization
	float max_rw=0.0, min_rw=1000.0, tot_rw=0.0;  //total rw
	float max_pc=0.0, min_pc=1000000.0, tot_pc=0.0;  //total pc
	float tmp_time=0;
	
	for (r=1; r<=nrun; r++) {
		int i, j, idle=0, unaloc_vm=0, ac_pm=0;
		float sum_cu=0, sum_ru=0, rw=0, pc=0, ave, var=0;
		struct machine pm_org[npm], pm_copy[npm], vm_copy[nvm];
		machines(npm, nvm);
		for (i=0; i<npm; i++) 
			pm_org[i] = pm_copy[i]= pm[i];
		for (i=0; i<nvm; i++) 
			vm_copy[i] = vm[i];
		
		gettimeofday(&start, NULL);
		
		vm_decreasing_order(vm_copy, nvm);	
		
		int min = 100000, rem, index;
		for (i=0; i<nvm; i++) {
			min = 100000;
			for (j=0; j<npm; j++) 
				if (vm_copy[i].c <= pm_org[j].c && vm_copy[i].r <= pm_org[j].r) {
					rem = pm_org[j].c - vm_copy[i].c;
					if (rem <min) {
						index = j;
						min = rem;
					}
				}
			pm_org[index].c -= vm_copy[i].c;
			pm_org[index].r -= vm_copy[i].r;
			vm_copy[i].c = 0;
		}
		
		gettimeofday(&end, NULL);
	    //printf("\nElapsed time: %ld microsecons\n", (end.tv_usec - start.tv_usec));
	    tmp_time += (end.tv_usec - start.tv_usec);
		
		for (i=0; i<npm; i++) {
			pm_org[i].cu = (float)(pm_copy[i].c-pm_org[i].c)/pm_copy[i].c;
			pm_org[i].ru = (float)(pm_copy[i].r-pm_org[i].r)/pm_copy[i].r;
			if (pm_org[i].c == pm_copy[i].c)
				idle++;
			else
				rw += (fabs(pm_org[i].cu-pm_org[i].ru)+0.0001)/(pm_org[i].cu+pm_org[i].ru);
		}
		for (i=0; i<npm; i++)
			if (pm_org[i].cu>0) {
				sum_cu += pm_org[i].cu;
				sum_ru += pm_org[i].ru;
				pc += (pm_org[i].p_min+(pm_org[i].p_max-pm_org[i].p_min)*pm_org[i].cu);
				ave = (pm_org[i].cu + pm_org[i].ru)/2;
				var += 0.5*((pm_org[i].cu - ave,2)+pow(pm_org[i].ru - ave,2));
			}
		for (i=0; i<nvm; i++) 
			if (vm_copy[i].c != 0)
				unaloc_vm++;
		ac_pm = npm - idle;
		
		// number of active pms	
		sum_ac_pm += ac_pm;
		if (ac_pm > max_ac_pm)
			max_ac_pm = ac_pm;
		if (ac_pm < min_ac_pm)
			min_ac_pm = ac_pm;
		// resource wastage
		tot_rw += rw;
		if (rw > max_rw)
			max_rw = rw;
		if (rw < min_rw)
			min_rw = rw;
		// power consumption 
		tot_pc += pc;
		if (pc > max_pc)
			max_pc = pc;
		if (pc < min_pc)
			min_pc = pc;
	}
	
	printf("\nResults:\n");
	printf("Algorithm: BFD, Run=%d\n", nrun);
	printf("Average number of PMs used: (max=%.2f , min=%.2f , ave=%.2f)\n",max_ac_pm-(float)sum_ac_pm/nrun, (float)sum_ac_pm/nrun-min_ac_pm, (float)sum_ac_pm/nrun); 
	printf("Average Resource Wastage: (max=%.2f , min=%.2f , ave=%.2f)\n",max_rw-tot_rw/nrun, tot_rw/nrun-min_rw, tot_rw/nrun);
	printf("Average Power Consumption: (max=%.2f , min=%.2f , ave=%.2f)\n\n",max_pc-tot_pc/nrun, tot_pc/nrun-min_pc, tot_pc/nrun);
}

//////////////////////////////  GRVMP  //////////////////////////////
void GRVMP(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int N=4;  // The number of Choices
	int max_unaloc_vm=0, min_unaloc_vm=10000, sum_unaloc_vm=0; // number of unallocated vms
	int max_ac_pm=0, min_ac_pm=10000, sum_ac_pm=0;  //number of active pms
	float max_cu=0.0, min_cu=1.0, tot_cu=0.0;  //total cpu utilization
	float max_ru=0.0, min_ru=1.0, tot_ru=0.0;  //total ram utilization
	float max_rw=0.0, min_rw=1000.0, tot_rw=0.0;  //total rw
	float max_pc=0.0, min_pc=1000000.0, tot_pc=0.0;  //total pc
	float tmp_time=0;
	
	for (r=1; r<=nrun; r++) {
		int i, j, idle=0, unaloc_vm=0, nvm1=nvm, x, flag, ac_pm=0;
		float sum_cu=0, sum_ru=0, rw=0, pc=0, ave, var=0, y;
		struct machine pm_org[npm], pm_copy[npm], vm_copy[nvm], vm1[nvm];
		machines(npm, nvm);
		
		for (i=0; i<npm; i++) 
			pm_org[i] = pm[i];
		for (i=0; i<nvm; i++) 
			vm_copy[i] = vm[i];			
		
		gettimeofday(&start, NULL);
		
		pm_decreasing_order_powerEfficiency(pm_org, npm);  //sorting PMs by their power effciency in decreasing order
		for (i=0; i<npm; i++) 
			pm_copy[i] = pm_org[i];
		
		float a[N][4], av, thr=0.0;
		while (nvm1>0) {
			a[0][0]= -1;
			for (i=0; i<(N<nvm1?N:nvm1); i++) {
				flag =1;
				x=rand()%nvm1;
				for(j=0; j<i; j++)
					if (x==a[j][0]) {
						i--;
						flag =0;
						break;
					}
				if (flag==0)
					continue;
				a[i][0]=x;
				a[i][1]=0; //cu
				a[i][2]=0; //ru
			}
			
			for (j=0; j<npm; j++) {
				for (i=0; i<(N<nvm1?N:nvm1); i++) {  //min(N,nvm1)
					x = (int)a[i][0];
					if (vm_copy[x].c <= pm_org[j].c && vm_copy[x].r <= pm_org[j].r) {
						a[i][1] = (float)(pm_copy[j].c-pm_org[j].c+vm_copy[x].c)/pm_copy[j].c;
						a[i][2] = (float)(pm_copy[j].r-pm_org[j].r+vm_copy[x].r)/pm_copy[j].r;
					}
				}
				flag = 0;
				for (i=0; i<(N<nvm1?N:nvm1); i++)
					if (a[i][1]!=0) {
						flag = 1;
						break;
					}	
				if (flag==0)
					continue;
				else {
					float min_rw=1000, max_ut=0;
					int min_rw_index,max_ut_index;
					for (i=0; i<(N<nvm1?N:nvm1); i++)
						if (a[i][1]!=0) {
							a[i][3] = (fabs(a[i][1]-a[i][2])+0.0001)/(a[i][1]+a[i][2]);	
							if (a[i][3]<min_rw) {
								min_rw=a[i][3];
								min_rw_index=i;
							}
							if (a[i][1]<thr && a[i][2]<thr && a[i][1]*a[i][2] > max_ut) {
								max_ut=a[i][1]*a[i][2];
								max_ut_index=i;
							}
						}
						if(max_ut > 0) {
							int index=a[max_ut_index][0];
							pm_org[j].c -= vm_copy[index].c;
							pm_org[j].r -= vm_copy[index].r;
							vm1[vm_copy[index].id].c=0;	
							for (int k=index; k<nvm1-1; k++)
								vm_copy[k]=vm_copy[k+1];
							nvm1--;
						}
						else {
							int index = a[min_rw_index][0];
							pm_org[j].c -= vm_copy[index].c;
							pm_org[j].r -= vm_copy[index].r;
							vm1[vm_copy[index].id].c=0;	
							for (int k=index; k<nvm1-1; k++)
								vm_copy[k]=vm_copy[k+1];
							nvm1--;
						}
					break;		
				}
			}
		}

		gettimeofday(&end, NULL);
	    //printf("\nElapsed time: %ld microsecons\n", (end.tv_usec - start.tv_usec));
		tmp_time += (end.tv_usec - start.tv_usec);
		
		for (i=0; i<npm; i++) {
			pm_org[i].cu = (float)(pm_copy[i].c-pm_org[i].c)/pm_copy[i].c;
			pm_org[i].ru = (float)(pm_copy[i].r-pm_org[i].r)/pm_copy[i].r;
			if (pm_org[i].c == pm_copy[i].c)
				idle++;
			else
				rw += (fabs(pm_org[i].cu-pm_org[i].ru)+0.0001)/(pm_org[i].cu+pm_org[i].ru);
		}
		
		for (i=0; i<npm; i++)
			if (pm_org[i].cu>0) {
				sum_cu += pm_org[i].cu;
				sum_ru += pm_org[i].ru;
				pc += (pm_org[i].p_min+(pm_org[i].p_max-pm_org[i].p_min)*pm_org[i].cu);
				ave = (pm_org[i].cu+pm_org[i].ru)/2;
				var += 0.5*((pm_org[i].cu - ave,2)+pow(pm_org[i].ru - ave,2));
			}
		for (i=0; i<nvm; i++) 
			if (vm1[i].c != 0)
				unaloc_vm++;
		ac_pm = npm - idle;
		
		// number of active pms	
		sum_ac_pm += ac_pm;
		if (ac_pm > max_ac_pm)
			max_ac_pm = ac_pm;
		if (ac_pm < min_ac_pm)
			min_ac_pm = ac_pm;
		// resource wastage
		tot_rw += rw;
		if (rw > max_rw)
			max_rw = rw;
		if (rw < min_rw)
			min_rw = rw;
		// power consumption 
		tot_pc += pc;
		if (pc > max_pc)
			max_pc = pc;
		if (pc < min_pc)
			min_pc = pc;
	}
	
	printf("\nResults:\n");
	printf("Algorithm: GRVMP, d=%d, Run=%d\n\n",N, nrun);
	printf("Average number of PMs used: (max=%.2f , min=%.2f , ave=%.2f)\n",max_ac_pm-(float)sum_ac_pm/nrun, (float)sum_ac_pm/nrun-min_ac_pm, (float)sum_ac_pm/nrun); 
	printf("Average Resource Wastage: (max=%.2f , min=%.2f , ave=%.2f)\n",max_rw-tot_rw/nrun, tot_rw/nrun-min_rw, tot_rw/nrun);
	printf("Average Power Consumption: (max=%.2f , min=%.2f , ave=%.2f)\n\n",max_pc-tot_pc/nrun, tot_pc/nrun-min_pc, tot_pc/nrun);
}	

//////////////////////////////  CRBFF  //////////////////////////////
void CRBFF(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int N=0;
	float tmp=(float)nvm/(float)npm;
//	printf("\n%.3f\n", tmp);
	if(ceil(tmp)<2)
		N=2;
	else
//		N=ceil(tmp)+round(nvm/npm);
		N=ceil(tmp);
//		N*=round(nvm/npm);
	
	int max_unaloc_vm=0, min_unaloc_vm=10000, sum_unaloc_vm=0; // number of unallocated vms
	int max_ac_pm=0, min_ac_pm=10000, sum_ac_pm=0;  //number of active pms
	float max_cu=0.0, min_cu=1.0, tot_cu=0.0;  //total cpu utilization
	float max_ru=0.0, min_ru=1.0, tot_ru=0.0;  //total ram utilization
	float max_rw=0.0, min_rw=1000.0, tot_rw=0.0;  //total rw
	float max_pc=0.0, min_pc=1000000.0, tot_pc=0.0;  //total pc
	float tmp_time=0;
	
	for (r=1; r<=nrun; r++) {
		int i, j, idle=0, unaloc_vm=0, nvm1=nvm, x, flag, ac_pm=0;
		float sum_cu=0, sum_ru=0, rw=0, pc=0, ave, var=0, y;
		struct machine pm_org[npm], pm_copy[npm], vm_copy[nvm], vm1[nvm];
		machines(npm, nvm);

		for (i=0; i<npm; i++) 
			pm_org[i] = pm[i];
		for (i=0; i<nvm; i++) 
			vm_copy[i] = vm[i];
			
		gettimeofday(&start, NULL);		
		
		pm_decreasing_order_powerEfficiency(pm_org, npm);  //sorting PMs by their power effciency in decreasing order
		for (i=0; i<npm; i++) 
			pm_copy[i] = pm_org[i];
		
		float a[N][4], av, thr=0.0;
		while (nvm1>0) {
			a[0][0]= -1;
			for (i=0; i<(N<nvm1?N:nvm1); i++) {
				flag =1;
				x=rand()%nvm1;
				for(j=0; j<i; j++)
					if (x==a[j][0]) {
						i--;
						flag =0;
						break;
					}
				if (flag==0)
					continue;
				a[i][0]=x;
				a[i][1]=0; //cu
				a[i][2]=0; //ru
			}
			
			for (j=0; j<npm; j++) {
				for (i=0; i<(N<nvm1?N:nvm1); i++) {  //min(N,nvm1)
					x = (int)a[i][0];
					if (vm_copy[x].c <= pm_org[j].c && vm_copy[x].r <= pm_org[j].r) {
						a[i][1] = (float)(pm_copy[j].c-pm_org[j].c+vm_copy[x].c)/pm_copy[j].c;
						a[i][2] = (float)(pm_copy[j].r-pm_org[j].r+vm_copy[x].r)/pm_copy[j].r;
					}
				}
				flag = 0;
				for (i=0; i<(N<nvm1?N:nvm1); i++)
					if (a[i][1]!=0) {
						flag = 1;
						break;
					}	
				if (flag==0)
					continue;
				else {
					float min_rw=1000, max_ut=0;
					int min_rw_index,max_ut_index;
					for (i=0; i<(N<nvm1?N:nvm1); i++)
						if (a[i][1]!=0) {
							a[i][3] = (fabs(a[i][1]-a[i][2])+0.0001)/(a[i][1]+a[i][2]);	
							if (a[i][3]<min_rw) {
								min_rw=a[i][3];
								min_rw_index=i;
							}
							if (a[i][1]<thr && a[i][2]<thr && a[i][1]*a[i][2] > max_ut) {
								max_ut=a[i][1]*a[i][2];
								max_ut_index=i;
							}
						}
						if(max_ut > 0) {
							int index=a[max_ut_index][0];
							pm_org[j].c -= vm_copy[index].c;
							pm_org[j].r -= vm_copy[index].r;
							vm1[vm_copy[index].id].c=0;	
							for (int k=index; k<nvm1-1; k++)
								vm_copy[k]=vm_copy[k+1];
							nvm1--;
						}
						else {
							int index = a[min_rw_index][0];
							pm_org[j].c -= vm_copy[index].c;
							pm_org[j].r -= vm_copy[index].r;
							vm1[vm_copy[index].id].c=0;	
							for (int k=index; k<nvm1-1; k++)
								vm_copy[k]=vm_copy[k+1];
							nvm1--;
						}
					break;		
				}
			}
		}
		
		gettimeofday(&end, NULL);
	    //printf("\nElapsed time: %ld microsecons\n", (end.tv_usec - start.tv_usec));
		tmp_time += (end.tv_usec - start.tv_usec);
		
		for (i=0; i<npm; i++) {
			pm_org[i].cu = (float)(pm_copy[i].c-pm_org[i].c)/pm_copy[i].c;
			pm_org[i].ru = (float)(pm_copy[i].r-pm_org[i].r)/pm_copy[i].r;
			if (pm_org[i].c == pm_copy[i].c)
				idle++;
			else
				rw += (fabs(pm_org[i].cu-pm_org[i].ru)+0.0001)/(pm_org[i].cu+pm_org[i].ru);
		}
		
		for (i=0; i<npm; i++)
			if (pm_org[i].cu>0) {
				sum_cu += pm_org[i].cu;
				sum_ru += pm_org[i].ru;
				pc += (pm_org[i].p_min+(pm_org[i].p_max-pm_org[i].p_min)*pm_org[i].cu);
				ave = (pm_org[i].cu+pm_org[i].ru)/2;
				var += 0.5*((pm_org[i].cu - ave,2)+pow(pm_org[i].ru - ave,2));
			}
		for (i=0; i<nvm; i++) 
			if (vm1[i].c != 0)
				unaloc_vm++;
		ac_pm = npm - idle;
		
		// number of active pms	
		sum_ac_pm += ac_pm;
		if (ac_pm > max_ac_pm)
			max_ac_pm = ac_pm;
		if (ac_pm < min_ac_pm)
			min_ac_pm = ac_pm;
		// resource wastage
		tot_rw += rw;
		if (rw > max_rw)
			max_rw = rw;
		if (rw < min_rw)
			min_rw = rw;
		// power consumption 
		tot_pc += pc;
		if (pc > max_pc)
			max_pc = pc;
		if (pc < min_pc)
			min_pc = pc;
	}
	
	printf("\nResults:\n");
	printf("Algorithm: CRBFF, N=%d, Run=%d\n",N, nrun);
	printf("Average number of PMs used: (max=%.2f , min=%.2f , ave=%.2f)\n",max_ac_pm-(float)sum_ac_pm/nrun, (float)sum_ac_pm/nrun-min_ac_pm, (float)sum_ac_pm/nrun); 
	printf("Average Resource Wastage: (max=%.2f , min=%.2f , ave=%.2f)\n",max_rw-tot_rw/nrun, tot_rw/nrun-min_rw, tot_rw/nrun);
	printf("Average Power Consumption: (max=%.2f , min=%.2f , ave=%.2f)\n\n",max_pc-tot_pc/nrun, tot_pc/nrun-min_pc, tot_pc/nrun);
}

//////////////////////////////  Functions  //////////////////////////////
void vm_decreasing_order(struct machine vm[], int nvm)
{
	int i , j;
	struct machine temp;
	for (i=0; i<nvm-1; i++)
		for (j=i+1; j<nvm; j++)
			if (vm[i].c < vm[j].c)
			{
				temp = 	vm[i];
				vm[i] = vm[j];
				vm[j] = temp;
			}
}

void pm_decreasing_order_powerEfficiency(struct machine pm[], int npm)
{
	int i, j;
	struct machine temp;
	for (i=0; i<npm; i++)
		pm[i].pow_eff = pm[i].c / pm[i].p_max;
		
	for (i=0; i<npm-1; i++)
		for (j=i+1; j<npm; j++) 
			if (pm[i].pow_eff < pm[j].pow_eff) {
				temp = 	pm[i];
				pm[i] = pm[j];
				pm[j] = temp;
			}
}

void machines(int npm, int nvm) {
	
	for (i=0; i<npm; i++) {
//================= GCP =================
//		pm[i].c = rand()%6001 + 2000;
//		pm[i].r = rand()%24577 + 8192;
//		pm[i].p_max = rand()%501 + 100;
//		pm[i].p_min = 0.7*pm[i].p_max;

//================= AWS EC2 =================
		pm[i].c = rand()%6901 + 4000;
		pm[i].r = rand()%10753 + 8192;
		pm[i].p_max = rand()%501 + 100;
		pm[i].p_min = 0.7*pm[i].p_max;

//================= Azure =================
//		pm[i].c = rand()%13601 + 6400;
//		pm[i].r = rand()%13313 + 7168;
//		pm[i].p_max = rand()%501 + 100;
//		pm[i].p_min = 0.7*pm[i].p_max;

		pm[i].id = i;
	}
//============================================

	for (i=0; i<nvm; i++) {	
//=============== GCP ===============
//		x = rand()%4;
//		if (x==0) {
//			vm[i].c = 250;
//			vm[i].r = 1024;
//		}	
//		else if (x==1) {
//			vm[i].c = 500;
//			vm[i].r = 2048;
//		}
//		else if (x==2) {
//			vm[i].c = 1000;
//			vm[i].r = 4096;
//		}
//		else  {
//			vm[i].c = 2000;
//			vm[i].r = 8192;
//		}

//=============== AWS EC2 ===============
		x = rand()%4;
		if (x==0) {
			vm[i].c = 500;
			vm[i].r = 613;
		}	
		else if (x==1) {
			vm[i].c = 1000;
			vm[i].r = 1740;
		}
		else if (x==2) {
			vm[i].c = 2000;
			vm[i].r = 3840;
		}
		else  {
			vm[i].c = 2500;
			vm[i].r = 870;
		}

//=============== Azure ===============
//		x = rand()%4;
//		if (x==0) {
//			vm[i].c = 1000;
//			vm[i].r = 768;
//		}	
//		else if (x==1) {
//			vm[i].c = 1600;
//			vm[i].r = 1792;
//		}
//		else if (x==2) {
//			vm[i].c = 3200;
//			vm[i].r = 3584;
//		}
//		else  {
//			vm[i].c = 6400;
//			vm[i].r = 7168;
//		}

		vm[i].id = i;
	}
}
