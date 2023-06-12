#ifndef	wynn_included
#define	wynn_included
#define	EIGEN_STACK_ALLOCATION_LIMIT	0
#define	EIGEN_NO_DEBUG
#include	<Eigen/Eigen>
#include	<cstdint>
#include	<cstdlib>
#include	<cstdio>
#include	<cfloat>
#include	<cmath>
#include	<ctime>
#include	<omp.h>
using	namespace	std;
using	namespace	Eigen;
namespace	wynn{
struct	global{
	unsigned	threads;	//	number of threads
	uint64_t	prng;	//	prng seed
	float	eta;	//	learning rate
	global(){	set_threads(omp_get_num_procs()>1?omp_get_num_procs()-1:1);	prng=time(NULL);	eta=1.0f/(1ull<<16);	}
	void	set_threads(unsigned	n){	threads=n;	omp_set_num_threads(threads);	setNbThreads(threads);	}
}global;
static	inline	uint64_t	wyrand(uint64_t	*seed){	*seed+=0xa0761d6478bd642full;	uint64_t	see1=*seed^0xe7037ed1a0b428dbull;	see1*=(see1>>32)|(see1<<32);	return	(*seed*((*seed>>32)|(*seed<<32)))^((see1>>32)|(see1<<32));	}
static	inline	double	wy2u01(uint64_t	r){	const	double	_wynorm=1.0/(1ull<<52);	return	(r>>12)*_wynorm;}
static	inline	float	wy2gau(uint64_t	r){	const	float	_wynorm=1.0/(1ull<<20);	return	((r&0x1fffff)+((r>>21)&0x1fffff)+((r>>42)&0x1fffff))*_wynorm-3.0f;}
static	inline	void	wymum(uint64_t	*A,	uint64_t	*B){	uint64_t	hh=(*A>>32)*(*B>>32),	hl=(*A>>32)*(uint32_t)*B,	lh=(uint32_t)*A*(*B>>32),	ll=(uint64_t)(uint32_t)*A*(uint32_t)*B;	*A=((hl>>32)|(hl<<32))^hh;	*B=((lh>>32)|(lh<<32))^ll;}
static	inline	uint64_t	wyhash64(uint64_t	A,	uint64_t	B){	A^=0xa0761d6478bd642full;	B^=0xa0761d6478bd642full;	wymum(&A,&B);	A^=0xa0761d6478bd642full;	B^=0xa0761d6478bd642full;	wymum(&A,&B);	return	A^B;}
template<unsigned	n>
static	inline	void	tiger(float	*w,	float	*g,	float	*m,	float	r){
	float	r2=r*2;
	#pragma	omp	parallel	for
	for(unsigned	i=0;	i<n;	i++){	m[i]+=0.03125f*(g[i]-m[i]);	w[i]-=r2*((m[i]>0)-0.5f);	}
}
}
#endif
