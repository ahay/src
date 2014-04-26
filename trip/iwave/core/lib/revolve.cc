/*----------------------------------------------------------------------------
 REVOLVE -- Checkpointing approaches
 File:     revolve.cpp
  
 Copyright (c) Andreas Griewank, Andrea Walther, Philipp Stumm
  
 This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.
----------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <list>
#include <algorithm>

using namespace std;

#include "revolve.h"

/************************************************************************************************************************************
All routines of class Schedule
**********************************************************************************************************************************/


Schedule::Schedule(int sn,Checkpoint *c)
{
	snaps=sn;
	checkpoint=c;
	checkpoint->ch[0]=0;
}




int Schedule::numforw(int steps, int snaps, ostream & str)
{

	int reps, range, num;

	if (snaps < 1)
	{
		str << " error occurs in numforw: snaps < 1 " << endl ;
		return -1;
	}
	if (snaps > checkup)
	{
		str << " number of snaps = " << snaps << " exceeds checkup " << endl;
		str << " redefine 'checkup'" << endl;
		return -1;
	}
	reps = 0;
	range = 1;
	while (range < steps)
	{
		reps += 1;
		range = range*(reps + snaps)/reps;
	}
	str << " range = " << range << " reps= " <<  reps << endl;
	if (reps > repsup)
	{
		str << " number of reps = " << reps << " exceeds repsup " << endl;
		str << " redefine 'repsup' " << endl;
		return -1;
	}
	num = reps * steps - range*reps/(snaps+1);
	return num;


}


/*************************************************************************************************************************************
All routines of class Online
***************************************************************************************************************************************/


Online::Online(int sn,Checkpoint *c,bool o=false) : Schedule(sn,c)
{

	checkpoint->init_ord_ch(); //  ord_ch = new int[snaps];
	//cout <<"\n Vectors size = " << checkpoint->ord_ch.size() << endl;
	//for(int i=0;i<snaps;i++) cout << " " << checkpoint->ord_ch[i];
	output=o;
	
}

Online::Online(Online &o) : Schedule(o.get_snaps())
{
	//ord_ch = new int[snaps];
	checkpoint = o.get_CP();
	//checkpoint->init_ord_ch();
	output=false;
	//advances = o.get_advances();
	//takeshots = o.get_shots();
	//commands = o.get_commands();
	check = o.get_check();
	capo = o.get_capo();

}

Online::~Online()
{
	//delete [] ord_ch;
	//delete [] num_rep;
}

/*************************************************************************************************************************************
All routines of class Online_r2
**************************************************************************************************************************************/

Online_r2::Online_r2(int sn,Checkpoint *c,bool o=false) : Online(sn,c,o)
{
	num_rep.reserve(snaps);
	check=-1;
	capo=0;
}

ACTION::action Online_r2::revolve(ostream & str)
{
	//cout <<" check = "<< check << " ch[check] " << checkpoint->ch[check]<< " capo " << capo << endl; 
	checkpoint->commands++;
	if ((check == -1) || ((checkpoint->ch[check] != capo) && (capo <= snaps-1)))
		// condition for takeshot for r=1
	{

		oldcapo_o = capo;
		check += 1;
		checkpoint->ch[check] = capo;
		t = 0;
		if (snaps < 4)
		{
			for(int i=0;i<snaps;i++)
				num_rep[i] = 2;
			incr = 2;
			iter = 1;
			oldind = snaps-1;
		}
		else
		{
			iter = 1;
			incr = 1;
			oldind = 1;
			for(int i=0;i<snaps;i++)
			{
				num_rep[i] = 1;
				checkpoint->ord_ch[i] = i;
			}
			offset = snaps-1;
		}
		if (capo == snaps-1)
		{
			ind = 2;
			old_f = 1;
		}
		/*   	  cout <<" 2 *check %d ch[*check] %d *capo %d endl",*check,ch[*check],*capo); */
		
		// Increase the number of takeshots and the corresponding checkpoint
		checkpoint->takeshots++;
		//checkpoint->number_of_writes[check]++;
		
		return ACTION::takeshot;
	}
	else if (capo < snaps-1)
		// condition for advance for r=1
	{
		capo = oldcapo_o+1;
		checkpoint->advances++;
		return ACTION::advance;
	}
	else
		//Online-Checkpointing for r=2
	{
		if (checkpoint->ch[check] == capo)
			// condition for advance for r=2
		{
			switch (snaps)
			{
				case 1: capo = MAXINT-1;
					checkpoint->advances++;
					return ACTION::advance;
				case 2: capo = checkpoint->ch[1]+incr;
					checkpoint->advances++;
					return ACTION::advance;
				case 3: checkpoint->advances+=incr;
					if (iter == 0)
					{
						capo = checkpoint->ch[oldind];
						for(int i=0;i<=(t+1)/2;i++)
						{
							capo += incr;
							incr++;
							iter++;
						}
					}
					else
					{
						capo = checkpoint->ch[ind]+incr;
						iter++;
						incr++;
					}
					
					return ACTION::advance;
				default: if (capo == snaps-1)
					{
						capo = capo+2;
						ind=snaps-1;
						checkpoint->advances+=2;
						return ACTION::advance;
					}
					if (output)
						str << " iter " << iter << " incr " << incr << "offset" << offset << endl;
					if (t == 0)
					{
						if (iter < offset)
						{
							capo = capo+1;
							checkpoint->advances++;
						}
						else
						{
							capo = capo+2;
							checkpoint->advances+=2;
						}
						if (offset == 1)
							t++;
						return ACTION::advance;
					}
					if (output)
						str << " iter " << iter << "incr " << incr << endl;
					str << " not implemented yet" << endl;
					return ACTION::error;
			}
		}
		else
			// takeshot for r=2
		{
			switch (snaps)
			{
				case 2: checkpoint->ch[1] = capo;
					incr++;
					// Increase the number of takeshots and the corresponding checkpoint
					checkpoint->takeshots++;
					//checkpoint->number_of_writes[1]++;
					return ACTION::takeshot;
				case 3: checkpoint->ch[ind] = capo;
					check = ind;
					str <<" iter " << iter << " num_rep[1] " << num_rep[1] << endl;
					if (iter == num_rep[1])
					{
						iter = 0;
						t++;
						oldind = ind;
						num_rep[1]++;
						ind = 2 - num_rep[1]%2;
						incr=1;
					}
					// Increase the number of takeshots and the corresponding checkpoint
					checkpoint->takeshots++;
					//checkpoint->number_of_writes[check]++;
					return ACTION::takeshot;
				default: if (capo < snaps+2)
					{
						checkpoint->ch[ind] = capo;
						check = ind;
						if (capo == snaps+1)
						{
							oldind = checkpoint->ord_ch[snaps-1];
							ind = checkpoint->ch[checkpoint->ord_ch[snaps-1]];
							if (output)
								str << " oldind " << oldind << " ind " << ind << endl;
							for(int k=snaps-1;k>1;k--)
							{
								checkpoint->ord_ch[k]=checkpoint->ord_ch[k-1];
								checkpoint->ch[checkpoint->ord_ch[k]] = checkpoint->ch[checkpoint->ord_ch[k-1]];
							}
							checkpoint->ord_ch[1] = oldind;
							checkpoint->ch[checkpoint->ord_ch[1]] = ind;
							incr=2;
							ind = 2;
							if (output)
							{
								str << " ind " << ind << " incr " << incr << " iter " << iter << endl;
								for(int j=0;j<snaps;j++)
									str << " j " << j << " ord_ch " << checkpoint->ord_ch[j] << " ch " << checkpoint->ch[checkpoint->ord_ch[j]] << " rep " << num_rep[checkpoint->ord_ch[j]] << endl;
							}
						}
						// Increase the number of takeshots and the corresponding checkpoint
						checkpoint->takeshots++;
						//checkpoint->number_of_writes[check]++;
						return ACTION::takeshot;
					}
					if (t == 0)
					{
						if (output)
							str <<" ind " << ind << " incr " <<  incr << " iter " << iter << " offset " << offset << endl;
						if (iter == offset)
						{
							offset--;
							iter = 1;
							check = checkpoint->ord_ch[snaps-1];
							checkpoint->ch[checkpoint->ord_ch[snaps-1]] = capo;
							oldind = checkpoint->ord_ch[snaps-1];
							ind = checkpoint->ch[checkpoint->ord_ch[snaps-1]];
							if (output)
								str << " oldind " << oldind << " ind " << ind << endl;
							for(int k=snaps-1;k>incr;k--)
							{
								checkpoint->ord_ch[k]=checkpoint->ord_ch[k-1];
								checkpoint->ch[checkpoint->ord_ch[k]] = checkpoint->ch[checkpoint->ord_ch[k-1]];
							}
							checkpoint->ord_ch[incr] = oldind;
							checkpoint->ch[checkpoint->ord_ch[incr]] = ind;
							incr++;
							ind=incr;
							if (output)
							{
								str << " ind " << ind << " incr " << incr << " iter " << iter << endl;
								for(int j=0;j<snaps;j++)
									str << " j " << j << " ord_ch " << checkpoint->ord_ch[j] << " ch " << checkpoint->ch[checkpoint->ord_ch[j]] << " rep " << num_rep[checkpoint->ord_ch[j]] <<  endl;
							}
						}
						else
						{
							checkpoint->ch[checkpoint->ord_ch[ind]] = capo;
							check = checkpoint->ord_ch[ind];
							iter++;
							ind++;
							if (output)
								str << " xx ind " << ind << " incr " << incr << " iter " << iter << endl;
						}
						// Increase the number of takeshots and the corresponding checkpoint
						checkpoint->takeshots++;
						//checkpoint->number_of_writes[check]++;
						return ACTION::takeshot;
					}
			}
		}
	}
	return ACTION::terminate;  // This means that the end of Online Checkpointing for r=2 is reached and
	//	another Online Checkpointing class must be started
}

/***********************************************************************************************************************************
All routines of class Online_r3
***********************************************************************************************************************************/

Online_r3::Online_r3(int sn,Checkpoint *c) : Online(sn,c)
{
	capo=(snaps+2)*(snaps+1)/2-1;
	ch3.reserve(snaps);
	cp_fest.reserve(snaps);
	tdiff.reserve(snaps);
	tdiff_end.reserve(snaps);
	check=1;
	for(int i=0;i<snaps;i++)
	{
		tdiff[i]=i+3;
		checkpoint->ord_ch[i]=snaps-i;
		cp_fest[i]=false;
	}
	tdiff_end[0]=6;
	for(int i=1;i<snaps;i++)
	{
		tdiff_end[i]=tdiff_end[i-1]+3+i;
	}
	ch3[0]=0;
	for(int i=1;i<snaps;i++)
	{
		ch3[i]=ch3[i-1]+tdiff_end[snaps-i-1];
	}

}

Online_r3::Online_r3(Online_r3 &o) : Online(o)
{
	capo=(snaps+2)*(snaps+1)/2-1;
	ch3.reserve(snaps);
	cp_fest.reserve(snaps);
	tdiff.reserve(snaps);
	tdiff_end.reserve(snaps);
	for(int i=0;i<snaps;i++)
	{
		tdiff[i]=i+3;
		//ord_ch[i]=snaps-i;
		cp_fest[i]=false;
	}
	tdiff_end[0]=6;
	for(int i=1;i<snaps;i++)
	{
		tdiff_end[i]=tdiff_end[i-1]+3+i;
	}
	ch3[0]=0;
	for(int i=1;i<snaps;i++)
	{
		ch3[i]=ch3[i-1]+tdiff_end[snaps-i-1];
	}
}


Online_r3::~Online_r3()
{
	//delete [] cp_fest;
	//delete [] tdiff;
	//delete [] tdiff_end;
	//delete [] ch3;
}

ACTION::action Online_r3::revolve(ostream & str)
{
	checkpoint->commands++;
	int n=1;
	if(capo==(snaps+2)*(snaps+1)/2-1)
		// Initialisation
	{

		capo+=1;
		forward=3;
		ind_now=1;
		checkpoint->advances+=3;
		cp=0;   /* changed 26.2.10*/
		return ACTION::advance;
	}
	else
	{
		if(capo==checkpoint->ch[check])
		{
			
			if(ind_now==snaps)
				forward = 1;
			else
			{
				if(capo==ch3[ind_now]-1)
					forward=1;
			}
			capo+=forward;
			checkpoint->advances+=forward;
			return ACTION::advance;
		}
		else if(capo<=(snaps+3)*(snaps+2)*(snaps+1)/6-4)
		{
			if(cp==0 && forward==1)
				//Now we take a checkpoint and the difference between the minimal number and this number
				// equals two
			{
				cp=0;
			}
			else
			{
				cp=choose_cp(n);
				while(cp_fest[checkpoint->ord_ch[snaps-1-cp]])
				{
					cp=choose_cp(++n);
				}
			}
			checkpoint->ch[checkpoint->ord_ch[snaps-1-cp]]=capo;
			tdiff_akt();
			akt_cp();
			check=checkpoint->ord_ch[snaps-1];
			if(checkpoint->ch[check]==ch3[ind_now])
				//saves a checkpoint that cannot be overwritten
			{
				cp_fest[check]=true;
				ind_now++;
			}
			forward=3;
			// Increase the number of takeshots and the corresponding checkpoint
			checkpoint->takeshots++;
			//checkpoint->number_of_writes[check]++;
			return ACTION::takeshot;
		}
		else
		{	// end of Online Checkpointing for r=3 is reached
			return ACTION::terminate;
		}
	}

}

// Function for selection of a checkpoint that can be replaced
// important for Online-Checkpointing for r=3
int Online_r3::choose_cp(int number)
{
	int i=2;
	if(tdiff[0]==3 && number==1) return 0;
	if(tdiff[0]+tdiff[1]<=10 && number<=2) return 1;
	while(number>0)
	{
		if(tdiff[i-1]+tdiff[i]<=tdiff_end[i])
		{
			number--;
		}
		i++;
	}
	return i-1;
}

// Renews the differences between the checkpoints
// Online-Checkpointing for r=3
void Online_r3::tdiff_akt()
{
	int i,sum;
	if(cp==0)
	{
		if(forward==3) tdiff[0]=6;
		else	tdiff[0]+=1;
		return;
	}
	else
	{
		sum=tdiff[0];
		//cp[0]=3;
	}
	for(i=cp-1;i>0;i--)
	{
		sum+=tdiff[i]-tdiff[i-1];
		tdiff[i]=tdiff[i-1];
	}
	tdiff[cp]+=sum;
	tdiff[0]=3;
}

// Renews the array of the indices of checkpoints
// Online-Checkpointing for r=3
void Online_r3::akt_cp()
{
	int i;
	if(cp==0) return;
	int value=checkpoint->ord_ch[snaps-1-cp];
	for(i=cp;i>0;i--)
	{
		//c[i+1]=c[i];
		checkpoint->ord_ch[snaps-i-1]=checkpoint->ord_ch[snaps-i];
	}
	//c[0]=value;
	checkpoint->ord_ch[snaps-1]=value;
}




/*******************************************************************************************************************************************
All routines of class Arevolve
********************************************************************************************************************************************/

Arevolve::Arevolve(int sn,Checkpoint *c, ostream & _str) : Online(sn,c), str(_str) 
{
	checkmax=snaps-1;
	capo=(snaps+3)*(snaps+2)*(snaps+1)/6-1;
	fine = capo+2;
	check=snaps-1;
	oldcapo=capo;
	newcapo=capo;
}

Arevolve::Arevolve(Arevolve &o) : Online(o), str(o.str)
{
	checkmax=snaps-1;
	fine = o.get_capo()+2;
	oldcapo=capo;
	newcapo=capo;
	check=snaps-1;
}


int Arevolve::tmin(int steps, int snaps, ostream & str)
{
	int reps, range, num;

	if (snaps < 1)
	{
		str << " error occurs in tmin: snaps < 1 " << endl;
		return -1;
	}
	if (snaps > checkup)
	{
		str << " number of snaps = " << snaps << " exceeds checkup " << endl;
		str << " redefine 'checkup' " << endl;
		return -1;
	}
	reps = 0;
	range = 1;
	while (range < steps)
	{
		reps += 1;
		range = range*(reps + snaps)/reps;
	}
	if (reps > repsup)
	{
		str << " number of reps = " << reps << " exceeds repsup " << endl;
		str << " redefine 'repsup' " << endl;
		return -1;
	}
	num = reps * steps - range*reps/(snaps+1);
	return num;
}

/* ************************************************************************* */

int Arevolve::sumtmin()
{
	int p=0, i;
	//if ( (check < 1)  &&  (steps>1) )
	//    return  MAXINT;
	for ( i=0; i<snaps-1; i++ )
	  p += tmin(checkpoint->ch[checkpoint->ord_ch[i+1]]-checkpoint->ch[checkpoint->ord_ch[i]],snaps-i,str);
	p = p + tmin(fine-1-checkpoint->ch[checkpoint->ord_ch[snaps-1]], 1, str) + fine-1;
	return p;
}

/* ************************************************************************* */

int Arevolve::mintmin( )
{
	int G=MAXINT, k=0, z=0, sum=0, g=0,i;


	sum=sumtmin();
	for (int j=1; j<snaps; j++)
	{
		g = z;
		if ( j-2>=0 )
		{
		  g = z+tmin( checkpoint->ch[checkpoint->ord_ch[j-1]]-checkpoint->ch[checkpoint->ord_ch[j-2]], snaps-j+2, str ) ;
			z=g;
		}
		if(j<snaps-1)
		{
		  g += tmin( checkpoint->ch[checkpoint->ord_ch[j+1]]-checkpoint->ch[checkpoint->ord_ch[j-1]], snaps-j+1, str ) ;
			for (i=j+1 ; i<=snaps-2; i++)
			  g += tmin( checkpoint->ch[checkpoint->ord_ch[i+1]]-checkpoint->ch[checkpoint->ord_ch[i]], snaps-i+1 , str) ;
			g+=tmin( fine-1-checkpoint->ch[checkpoint->ord_ch[snaps-1]], 2, str ) ;
		}
		else
		  g+=tmin( fine-1-checkpoint->ch[checkpoint->ord_ch[snaps-2]], 2 , str);
		if (g < G )
		{
			G = g;
			k = j;
		}
	}
	g = G + fine-1;
	if ( g < sum )
		return k;
	else
		return 0;
}

/* ************************************************************************* */
void Arevolve::akt_cp(int cp)
{
	int value =checkpoint->ord_ch[cp];

	for(int j=cp;j<snaps-1;j++)
		checkpoint->ord_ch[j]=checkpoint->ord_ch[j+1];
	checkpoint->ord_ch[snaps-1]=value;
}


ACTION::action Arevolve::revolve(ostream & str)
{

	oldcapo=capo;
	int shift = mintmin( );
	checkpoint->commands++;

	if ( shift==0 )
	{
		capo = oldcapo+1;
		oldfine = fine;
		//advances += capo-oldcapo;
		fine++;
		checkpoint->advances++;
		return  ACTION::advance;   //while arevolve
	}
	else
	{
		capo = oldcapo+1;
		checkpoint->ch[checkpoint->ord_ch[shift]]=capo;
		akt_cp(shift);
		check=checkpoint->ord_ch[shift];
		oldfine = fine++;
		newcapo = capo;
		checkpoint->takeshots++;
		//checkpoint->number_of_writes[check]++;
		return   ACTION::takeshot;
	}

}

/*******************************************************************************************************************************************
All routines of class Moin
********************************************************************************************************************************************/

Moin::Moin(int sn,Checkpoint *c) : Online(sn,c)
{
	capo=(snaps+3)*(snaps+2)*(snaps+1)/6-1;
	d.reserve(snaps);
	l.reserve(snaps);
	l[0]=10000;
	d[0]=false;
	for(int i=1;i<snaps;i++)
	{
		l[i]=2;
		d[i]=true;
	}
	start=true;
	start1=true;
	is_d=false;
	check=0;
}

Moin::Moin(Moin &o) : Online(o)
{
	
}

bool Moin::is_dispensable(int *index)
{
	bool dis=false;
	int ind=0;
	for(int i=snaps-1;i>0;i--)
	{
		if(d[i])
		{
			dis=true;
			if(checkpoint->ch[i]>ind)
			{
				ind=checkpoint->ch[i];
				*index = i;
			}
		}
	}
	return dis;
		
}

int Moin::get_lmin()
{
	int lmin=l[1];

	for(int i=2;i<snaps;i++)
	{
		if(l[i]<lmin)
			lmin=l[i];
	}
	return lmin;
}

void Moin::adjust_cp(int index)
{
	int level = l[index];
	int time  = checkpoint->ch[index];

	for(int i=snaps-1;i>0;i--)
	{
		if(i != index)
		{
			if(l[i]<level && checkpoint->ch[i] < time)
			{
				d[i] = true;
			}
		}
	}
}

ACTION::action Moin::revolve(ostream & str)
{
	int index;
	//checkpoint->print_ch(cout);
	//checkpoint->print_ord_ch(cout);
	
	//cout << "\n \n capo = " << capo;
	//cout << "\n \n snaps = " << snaps;
	checkpoint->commands++;
	if(start) 
	{
		capo++;
		start=false;
		checkpoint->advances++;
		return ACTION::advance;
	}
	if(start1)
	{
		start1=false;
		for(int i=1;i<snaps;i++)
		{
			if(checkpoint->ord_ch[i]==snaps-1)
			{
				checkpoint->ch[i] = capo;
				check=i;
				l[i]=3;
				d[i]=false;
				
			}
		}	
		forward=1;
		// Increase the number of takeshots and the corresponding checkpoint
		checkpoint->takeshots++;
		//checkpoint->number_of_writes[check]++;
		return ACTION::takeshot;
	}
	if(forward>0)
	{
		capo +=forward;
		forward=0;
		checkpoint->advances++;
		return ACTION::advance;
	}
	else
	{
		if(is_dispensable(&index))
		{
			//cout << "is_dispensable " << endl;
			checkpoint->ch[index] = capo;
			l[index] = 0;
			d[index] = false;
			index_old = index;	
			forward=1;
			check = index;
			checkpoint->takeshots++;
			return ACTION::takeshot;
		}
		else if(is_d)
		{
			//cout << "is_d=true " << endl;
			checkpoint->ch[index_old] = capo;
			check=index_old;
			lmin = get_lmin();
			l[index_old] = lmin+1;
			d[index_old]=false;
			//cout << "check = " << check << "  forward= " << lmin+1 << endl;
			adjust_cp(index_old);
			is_d=false;
			forward=1;
			checkpoint->takeshots++;
			return ACTION::takeshot;
		}
		else
		{
			lmin = get_lmin();
			//cout << "lmin = " << lmin << endl;
			forward = lmin+1;
			capo +=forward;
			is_d=true;
			forward=0;
			checkpoint->advances++;
			return ACTION::advance;
			
			//cout << "i = " << i << endl;
			
		}
	}
	str << "\n \n Irgendwas ist falsch \n\n";
	return ACTION::terminate;	
}

/***************************************************************************************************************************************
All routines of class Offline
****************************************************************************************************************************************/



Offline::Offline(int st,int sn,Checkpoint *c, ostream & _str) : Schedule(sn,c), str(_str)
{
	//advances = 0;
	//takeshots = 0;
	//commands = 0;
	steps=st;
	check=-1;
	info=3;
	fine=steps;
	capo=0;
	online=false;
}

Offline::Offline( int sn,Checkpoint *c,Online *o,int f, ostream & _str) : Schedule(sn,c ), str(_str)
{

	checkpoint->ch[0]=0;
	online=true;
	check=o->get_check();

	capo=o->get_capo();
	turn=0;
	num_ch.reserve(snaps);
	ind = 0;
	for(int i=0;i<snaps;i++)
	{
		num_ch[i] = 0;
		for(int j=0;j<snaps;j++)
		{
			if (checkpoint->ch[j] < checkpoint->ch[i])
				num_ch[i]++;
		}
		if (o->get_output())
			str << " i " << i << " num_ch " << num_ch[i] << " ch " << checkpoint->ch[i] << endl;
		/*for(int k=0;k<snaps;k++)
		{
			for(int j=0;j<snaps;j++)
				if (num_ch[j] == k)
					checkpoint->ord_ch[k] = j;  
			//if (o->get_output())
				printf(" i %d ord_ch %d ch %d\n",k,checkpoint->ord_ch[k],checkpoint->ch[k]);
		}*/
	}
	for(int i=0;i<snaps;i++)
	{
		for(int j=0;j<snaps;j++)
			if (num_ch[j] == i)
				checkpoint->ord_ch[i] = j; 
		if (o->get_output())
			printf(" i %d ord_ch %d ch %d\n",i,checkpoint->ord_ch[i],checkpoint->ch[i]);
	}
	checkpoint->advances = f-1;
	info=3;
	//takeshots = o->get_shots();
	//commands = o->get_commands();
	oldsnaps=snaps;
}

Offline::Offline( Schedule *o, ostream & _str) : Schedule(o->get_snaps() ), str(_str) 
{}

Offline::Offline( Offline &o) : Schedule(o.get_snaps(),o.get_CP()), str(o.str)
{
	//advances = o.get_advances();
	//takeshots = o.get_shots();
	//commands = o.get_commands();
	steps=o.get_steps();
	check=o.get_check();
	info=3;
	fine=o.get_steps();
	capo=o.get_capo();
	online=o.get_online();
	num_ch.reserve( snaps);
	for(int i=0;i<snaps;i++)
	{
		num_ch[i]=o.get_num_ch(i);
	}
	oldsnaps=snaps;

}

ACTION::action Offline::revolve(ostream & str)
{
	checkpoint->commands++;
	if ((check < -1) || (capo > fine))
		return ACTION::error;
	if ((check == -1) && (capo < fine))
	{
		if (check == -1)
			turn = 0;   /* initialization of turn counter */
		checkpoint->ch[0] = capo-1;
	}
	switch (fine-capo)
	{
		case 0:   /* reduce capo to previous checkpoint, unless done  */
			if (check == -1 || capo==checkpoint->ch[0] )
			{
				if (info > 0)
				{
					str << "\n advances: " << setw(5) << checkpoint->advances;
					str << "\n takeshots: " << setw(5) << checkpoint->takeshots;
					str << "\n commands: " << setw(5) << checkpoint->commands << endl;
				}
				return ACTION::terminate;
			}
			else
			{
				if(online)
				{
					int ind = 0;
					for(int i=0;i<snaps;i++)
					{
						if ((checkpoint->ch[i] > checkpoint->ch[ind]) && (checkpoint->ch[i] < capo))
							ind = i;
					}
					check = ind;
				}
				capo = checkpoint->ch[check];
				oldfine = fine;
				checkpoint->number_of_reads[check]++;
				return ACTION::restore;
			}
		case 1:  /* (possibly first) combined forward/reverse step */
			fine -= 1;
			if (check >= 0 && checkpoint->ch[check] == capo)
				check -= 1;
			if (turn == 0)
			{
				turn = 1;
				oldfine = fine;
				return ACTION::firsturn;
			}
			else
			{
				oldfine = fine;
				return ACTION::youturn;
			}
		default:
			if (check==-1)
				// Initialisation
			{
				checkpoint->ch[0]=0;
				check=0;
				oldsnaps = snaps;
				if (snaps > checkup)
				{
					info = 14;
					return ACTION::error;
				}
				if (info > 0)
				{
				  int num = numforw(fine-capo,snaps,str);
					if (num == -1)
					{
						info = 12;
						return ACTION::error;
					}
					str << " prediction of needed forward steps: " << setw(8) << num << " => " << endl;
					str << " slowdown factor: " << setiosflags(ios::fixed) << setprecision(4) << ((double) num)/(fine-capo)<< endl << endl;
				}
				oldfine = fine;
				//last_action=takeshot;
				checkpoint->number_of_writes[check]++;
				checkpoint->takeshots++;
				return ACTION::takeshot;
			}
			if (checkpoint->ch[check]!=capo)
				// takeshot
			{
				if (online)					
					check=checkpoint->ord_ch[num_ch[check]+1];
				else
					check++;
				if (check >= checkup)
				{
					info = 10;
					return ACTION::error;
				}
				if (check+1 > snaps)
				{
					info = 11;
					return ACTION::error;
				}
				checkpoint->ch[check] = capo;
				checkpoint->takeshots++;
				oldfine = fine;
				//last_action=takeshot;
				checkpoint->number_of_writes[check]++;
				return ACTION::takeshot;
			}
			else
				// advance
			{
				if ((oldfine < fine) && (snaps == check+1))
				{
					info = 13;
					return ACTION::error;
				}
				int oldcapo = capo;
				int ds;
				if (online)
					ds = snaps - num_ch[check];
				else
					ds = snaps - check;
				if (ds < 1)
				{
					info = 11;
					return ACTION::error;
				}
				int reps = 0;
				int range = 1;
				while (range < fine - capo)
				{
					reps += 1;
					range = range*(reps + ds)/reps;
				}
				if (reps > repsup)
				{
					info = 15;
					return ACTION::error;
				}
				if (snaps != oldsnaps)
				{
					if (snaps > checkup)
					{
						info = 14;
						return ACTION::error;
					}
				}
				int bino1 = range*reps/(ds+reps);
				int bino2 = (ds > 1) ? bino1*ds/(ds+reps-1) : 1;
				int bino3;
				if (ds == 1)
					bino3 = 0;
				else
					bino3 = (ds > 2) ? bino2*(ds-1)/(ds+reps-2) : 1;
				int bino4 = bino2*(reps-1)/ds;
				int bino5;
				if (ds < 3)
					bino5 = 0;
				else
					bino5 = (ds > 3) ? bino3*(ds-2)/reps : 1;

				if (fine-capo <= bino1 + bino3)
					capo = capo+bino4;
				else
				{
					if (fine-capo >= range - bino5)
						capo = capo + bino1;
					else
						capo = fine-bino2-bino3;
				}
				if (capo == oldcapo)
					capo = oldcapo+1;
				checkpoint->advances += capo - oldcapo;
				oldfine = fine;
				//last_action=advance;
				return ACTION::advance;
			}
	}
}


// All routines of class Revolve

Revolve::Revolve(int st,int sn, ostream & _str)
  : str(_str) {
	capo=0;
	checkpoint = new Checkpoint(sn);
	f=new Offline(st,sn,checkpoint,str);
	online=false;
	steps=st;
	snaps=sn;
	check=-1;
	info = 1;
	multi=false;
	where.reserve(snaps);
	for(int i=0;i<snaps;i++)
		where[i] = true;
	checkpoint->advances=0;
	checkpoint->takeshots=0;
	checkpoint->commands=0;
}

Revolve::Revolve(int st,int sn,int sn_ram, ostream & _str)
  : str(_str) {
	vector <int> v;
	int num=0,mid;
	
	checkpoint = new Checkpoint(sn);
	f=new Offline(st,sn,checkpoint,str);
	online=false;
	steps=st;
	snaps=sn;
	check=-1;
	info = 1;
	multi=true;
	where.reserve(snaps);
	indizes_ram.reserve(snaps);
	indizes_rom.reserve(snaps);	
	
	v = get_write_and_read_counts();
	sort(v.begin(),v.end());
	mid=v[snaps-sn_ram];
	//cout << mid << endl;
	for(int i=snaps-1;i>=0;i--)
	{
		if(v[i]>=mid && num<sn_ram)
		{
			where[i]=true;
			num++;
		}
		else
		{
			where[i]=false;
		}
	}
	int j=0,k=0;
	for(int i=0;i<snaps;i++)
	{
		if(where[i])
			indizes_ram[i]=j++;
		else
			indizes_rom[i]=k++;
	} 		
	checkpoint->advances=0;
	checkpoint->takeshots=0;
	checkpoint->commands=0;
}

Revolve::Revolve(int sn, ostream & _str)
  : str(_str) {
	checkpoint = new Checkpoint(sn);
	f=new Online_r2(sn,checkpoint);
	where.reserve(sn);
	online=true;
	snaps=sn;
	//info=inf;
	check=-1;
	info = 1;
	r=2;
	checkpoint->advances=0;
	checkpoint->takeshots=0;
	checkpoint->commands=0;
}

ACTION::action Revolve::revolve(int* check,int* capo,int* fine,int snaps,int* info,bool *where_to_put)
{
	ACTION::action whatodo;
	whatodo=f->revolve(str);
	if(online && whatodo==ACTION::terminate && r==2)
	{
		f = new Online_r3(snaps,checkpoint);
		whatodo=f->revolve(str);
		r++;
	}
	if(online && whatodo==ACTION::terminate && r==3)
	{
		f = new Moin(snaps,checkpoint);
		whatodo=f->revolve(str);
		r++;
	}
	*check=f->get_check();
	*capo=f->get_capo();
	*fine=f->get_fine();
	//*info=f->get_info();
	if(*check==-1) return whatodo;
	*where_to_put=where[*check];
	return whatodo;

}

ACTION::action Revolve::revolve(int* check,int* capo,int* fine,int snaps,int* info)
{
	ACTION::action whatodo;
	whatodo=f->revolve(str);
	if(online && whatodo==ACTION::terminate && r==2)
	{
		delete f;
		f = new Online_r3(snaps,checkpoint);
		whatodo=f->revolve(str);
		r++;
	}
	if(online && whatodo==ACTION::terminate && r==3)
	{
		delete f;
		f = new Moin(snaps,checkpoint);
		whatodo=f->revolve(str);
		r++;
	}
	*check=f->get_check();
	*capo=f->get_capo();
	*fine=f->get_fine();
	//*info=f->get_info();
	return whatodo;
}

ACTION::action Revolve::revolve()
{
	ACTION::action whatodo;
	oldcapo = capo;
	whatodo=f->revolve(str);
	if(online && whatodo==ACTION::terminate && r==2)
	{
		delete f;
		f = new Online_r3(snaps,checkpoint);
		whatodo=f->revolve(str);
		r++;
	}
	if(online && whatodo==ACTION::terminate && r==3)
	{
		delete f;
		f = new Moin(snaps,checkpoint);
		whatodo=f->revolve(str);
		r++;
	}
	check=f->get_check();
	capo=f->get_capo();
	fine=f->get_fine();
	info=f->get_info();
	if(check<=-1) return whatodo;
	if(!online) where_to_put=where[check];
	return whatodo;
}




void Revolve::turn(int final)
{
	if(online)
	{
		fine=final;
		capo=final-1;
		Online *g = new Online((Online &) *f);
		delete f;
		f=new Offline(snaps,checkpoint,g,final,str);
		delete g;
		//f=new Offline(snaps,checkpoint,(Online*)f,final);
		f->set_fine(final);
		f->set_capo(final-1);
		online=false;
	}
}

double Revolve::expense(int steps, int snaps, ostream & str)
{
	double ratio;

	if (snaps < 1)
	{
		str << " error occurs in expense: snaps < 0 " << endl;
		return -1;
	}
	if (steps < 1)
	{
		str <<" error occurs in expense: steps < 0 " << endl;
		return -1;
	}
	ratio = ((double) numforw(steps,snaps,str));
	if (ratio == -1)
		return -1;
	ratio = ratio/steps;
	return ratio;
}

/* ************************************************************************* */

int Revolve::maxrange(int ss, int tt, ostream & str)
{
	int i, ires;
	double res = 1.0;

	if ((tt<0) || (ss<0))
	{
		str << " error in MAXRANGE: negative parameter ";
		return -1;
	}
	for (i=1; i<= tt; i++)
	{
		res *= (ss + i);
		res /= i;
		if (res > MAXINT)
		{
			ires=MAXINT;
			str << " warning from MAXRANGE: returned maximal integer "<< ires << endl;
			return ires;
		}
	}
	ires = (int) res;
	return ires;
}

/* ************************************************************************* */

int Revolve::adjust(int steps, ostream & str)
{
	int snaps, s, reps;

	snaps = 1;
	reps = 1;
	s = 0;
	while ( maxrange(snaps+s, reps+s, str) > steps )
		s--;
	while ( maxrange(snaps+s, reps+s, str) < steps )
		s++;
	snaps += s;
	reps += s ;
	s = -1;
	while ( maxrange(snaps,reps,str) >= steps )
	{
		if (snaps > reps)
		{
			snaps -= 1;
			s = 0;
		}
		else
		{
			reps -= 1;
			s = 1;
		}
	}
	if ( s == 0 )
		snaps += 1 ;
	if ( s == 1 )
		reps += 1;
	return snaps;
}

int Revolve::get_r(int steps,int snaps, ostream & str)
{
  int reps, range;
  //, num;

	if (snaps < 1)
	{
		str << " error occurs in tmin: snaps < 1 " << endl;
		return -1;
	}
	if (snaps > checkup)
	{
		str << " number of snaps = " << snaps << " exceeds checkup " << endl;
		str << " redefine 'checkup' " << endl;
		return -1;
	}
	reps = 0;
	range = 1;
	while (range < steps)
	{
		reps += 1;
		range = range*(reps + snaps)/reps;
	}
	if (reps > repsup)
	{
		str << " number of reps = " << reps << " exceeds repsup " << endl;
		str << " redefine 'repsup' " << endl;
		return -1;
	}
	return reps;
}

int Revolve::get_r( ostream & str )
{
  int reps, range;
  //, num;

	if (snaps < 1)
	{
		str << " error occurs in tmin: snaps < 1 " << endl;
		return -1;
	}
	if (snaps > checkup)
	{
		str << " number of snaps = " << snaps << " exceeds checkup " << endl;
		str << " redefine 'checkup' " << endl;
		return -1;
	}
	reps = 0;
	range = 1;
	while (range < steps)
	{
		reps += 1;
		range = range*(reps + snaps)/reps;
	}
	if (reps > repsup)
	{
		str << " number of reps = " << reps << " exceeds repsup " << endl;
		str << " redefine 'repsup' " << endl;
		return -1;
	}
	return reps;
}



vector <int> Revolve::get_write_and_read_counts()
{
	vector <int> num(snaps);
	
	for(int i=0;i<snaps;i++)
		num[i]=get_number_of_writes_i( steps,snaps,i) + get_number_of_reads_i(steps,snaps,i);
	
	return num;
}

int cal(int l,int c,int i)
{
	if(i==0) return 0;
	if(l>(1+i)*c-0.5*(i-1)*i+1) return i;
	return (int) floor(0.5*(1+2*c)-sqrt(pow(0.5*(1+2*c),2)+2*i+4-2*l))-1;
}


/***************************************************************************************************************************
* 	Determination of the write counts for the checkpoint i with l steps and c checkpoints
*
*
*	The write counts are determined "Multi-Stage Approaches for Optimal Offline Checkpointing"
*
*	min l			max l				write counts 		
*	
*	0			1+i				0 (Lemma 3.3)
*	2+i			2c+i				1 (Lemma 3.3)
*	beta(c,r-1)		beta(c,r-1)+beta(c-1,r-1)	beta(i,r-2) (Theorem 3.2)
*	beta(c,1)		beta(c,1)+beta(c-1,1)		i+1 (Theorem 3.2)
*	beta(c,2)+beta(c-1,2)	beta(c,3)			Algorithm I
*
***************************************************************************************************************************/




int Revolve::get_number_of_writes_i(int l, int c,int i)
{
	if (i==0) return 1;
	if(l <= 1+i) return 0;
	else if(l <= 2*c+i) return 1;
	else if(l <= (1+i)*c-0.5*(i-1)*i+1) return (int) floor(0.5*(1+2*c)-sqrt(pow(0.5*(1+2*c),2)+2*i+4-2*l));
	else if(l <=c*c+2.*c+i) return i+1;
	else
	{
		double l_0=c*c+2.*c+1.;
		double a=27.*c*(c*c-1.)+162.*(l_0-l);
		int k;

		if (a==0)	
			k=c-1;
		else
			k = (int) floor(c-pow(2./(a+sqrt(a*a-108.)),1./3.)-1./3.*pow(0.5*(a+sqrt(a*a-108.)),1./3.));
		
		double l_k = 1./6.*k*k*k-c/2.*k*k+1./6.*(3.*c*c-1.)*k+l_0;
		
		if(i<=k)
			return (int) (1./2.*i*i+3./2.*i+1.);
		else    
		{
			int w_i_k= (int) (i*k+i+1.- 1./2.*k*(k-1.));
			return w_i_k+cal(l-(int)l_k+2*(c-k)+1,c-k,i-k);
		}
	}
}

	
/***************************************************************************************************************************
* 	Determination of the number of read operations for the checkpoint i with l steps and c checkpoints
*
*
*	The read counts are determined "Multi-Stage Approaches for Optimal Offline Checkpointing"
*
*	min l				max l			# reads 		
*	
*	0				1+i			0 (Lemma 4.2)
*	2+i				2c-i			1 (Lemma 4.2)
*	beta(c,r-1)+beta(c-1,r-1)	beta(c,r)		w_i(c,l)+beta(i+1,r-2)	Theorem 4.5
*	2c+1				beta(c,2)		w_i(c,l)+1 Theorem 4.5
*	beta(c,2)+1			beta(c,2)+beta(c-1,2)	Theorem 4.8
*	beta(c,2)+beta(c-1,2)+1 	beta(c,3)		Theorem 4.5
*
***************************************************************************************************************************/
	
	
int Revolve::get_number_of_reads_i(int l, int c,int i)
{
	if(l <= i+1) return 0;
	else if (l <= 2*c-i) return 1;
	else if (l <= 2*c+1) return 2;
	else if (l <= c*c/2.+3.*c/2.+1) return get_number_of_writes_i(l,c,i)+1;
	else if (l <= c*c+2*c+1)	
	{	
		double a = 0.5*(4.*i-2.*c+7.);
		double b = pow(c-2.*i-3.,2.)+c+3.;
		if(l >= c*c+2*c+1-c*i+0.5*(i*i-i))
			return (int) floor(a+sqrt(a*a-b+2.*(l-0.5*c*c-1.5*c-1)));
		else
			return i+2;
	}
	else return get_number_of_writes_i(l,c,i)+i+2; 
}
