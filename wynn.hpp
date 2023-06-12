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
using	namespace	std;
using	namespace	Eigen;
namespace	wynn{
uint64_t	prng=time(NULL);
static	inline	uint64_t	wyrand(uint64_t	*seed){	*seed+=0xa0761d6478bd642full;	uint64_t	see1=*seed^0xe7037ed1a0b428dbull;	see1*=(see1>>32)|(see1<<32);	return	(*seed*((*seed>>32)|(*seed<<32)))^((see1>>32)|(see1<<32));	}
static	inline	double	wy2u01(uint64_t	r){	const	double	_wynorm=1.0/(1ull<<52);	return	(r>>12)*_wynorm;}
static	inline	float	wy2gau(uint64_t	r){	const	float	_wynorm=1.0/(1ull<<20);	return	((r&0x1fffff)+((r>>21)&0x1fffff)+((r>>42)&0x1fffff))*_wynorm-3.0f;}
static	inline	void	wymum(uint64_t	*A,	uint64_t	*B){	uint64_t	hh=(*A>>32)*(*B>>32),	hl=(*A>>32)*(uint32_t)*B,	lh=(uint32_t)*A*(*B>>32),	ll=(uint64_t)(uint32_t)*A*(uint32_t)*B;	*A=((hl>>32)|(hl<<32))^hh;	*B=((lh>>32)|(lh<<32))^ll;}
static	inline	uint64_t	wyhash64(uint64_t	A,	uint64_t	B){	A^=0xa0761d6478bd642full;	B^=0xa0761d6478bd642full;	wymum(&A,&B);	A^=0xa0761d6478bd642full;	B^=0xa0761d6478bd642full;	wymum(&A,&B);	return	A^B;}


}
#endif
