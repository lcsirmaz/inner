/** round.h -- rounding a double using continued fraction **/

/***********************************************************************
 * This code is part of INNER, a linear multiobjective problem solver.
 *
 * Copyright (C) 2016-2024 Laszlo Csirmaz, https://github.com/lcsirmaz/inner
 *
 * This program is free, open-source software. You may redistribute it
 * and/or modify under the terms of the GNU General Public License (GPL).
 *
 * There is ABSOLUTELY NO WARRANTY, use at your own risk.
 ***********************************************************************/
         
/**********************************************************************
* The number is exapnded as a countinued fraction using four iterations.
* If the error is smaller that PARAMS(RoundEps), then it is replaced by
* the value of the fraction.
*
* static inline long intfloor(double d)
*  an auxiliary procedure computing round(d)
*
* static inline void round_to(double *v)
*  the rounding procedure
*/
#ifndef INNER_ROUND_H
#define INNER_ROUND_H
#define ROUND_EPS	PARAMS(RoundEps)

static inline long intfloor(double d)
{ return d<0.0 ? (long)(d-0.5) : (long)(d+0.5); }

static inline void round_to(double *v)
{double ip,iip,v2,ip2,ip3,ip4;
    ip=intfloor(*v); iip=*v-ip;
//*v  (-2.5,-1.5]  (-1.5,-0.5] (-0.5,0] [0,0.5) [0.5,1.5)  [1.5,2.5)
//ip    -2          -1           0          0      1          2
//iip (-0.5,0.5]   (-0.5,0.5]  (-0.5,0] [0.0.5] [-0.5,0.5) [-0.5,0.5)
    if(iip<=-0.5||iip>=0.5) return;
    if(-ROUND_EPS<iip && iip<ROUND_EPS){ *v=ip; return; }
    // *v==ip+iip, -0.50<=iip<=0.5
    v2=1.0/iip; ip2=intfloor(v2); iip=v2-ip2;
    if(-2.0*ROUND_EPS<iip && iip<2.0*ROUND_EPS){ *v=ip+1.0/ip2; return; }
    v2=1.0/iip; ip3=intfloor(v2); iip=v2-ip3;
    if(-4.0*ROUND_EPS<iip && iip<4.0*ROUND_EPS){ *v=ip+1.0/(ip2+1.0/ip3); return; }
    v2=1.0/iip; ip4=intfloor(v2); iip=v2-ip4;
    if(-8.0*ROUND_EPS<iip && iip<8.0*ROUND_EPS){ *v=ip+1.0/(ip2+1.0/(ip3+1.0/ip4)); return;}
}

#undef ROUND_EPS
#endif /* INNER_ROUND_H */

/* EOF */

