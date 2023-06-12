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

namespace	wynn{

struct	global{
	uint64_t	prng,threads;
	float	learning_rate;
	global(){	set_threads(omp_get_num_procs()>1?omp_get_num_procs()/2:1);	prng=time(NULL);	learning_rate=1.0f/(1ull<<16);	}
	void	set_threads(uint64_t	n){	threads=n;	omp_set_num_threads(threads);	Eigen::setNbThreads(threads);	}
}global;

static	inline	uint64_t	wyrand(uint64_t	*seed){	*seed+=0xa0761d6478bd642full;	uint64_t	see1=*seed^0xe7037ed1a0b428dbull;	see1*=(see1>>32)|(see1<<32);	return	(*seed*((*seed>>32)|(*seed<<32)))^((see1>>32)|(see1<<32));	}

static	inline	double	wy2u01(uint64_t	r){	const	double	_wynorm=1.0/(1ull<<52);	return	(r>>12)*_wynorm;}

static	inline	float	wy2gau(uint64_t	r){	const	float	_wynorm=1.0/(1ull<<20);	return	((r&0x1fffff)+((r>>21)&0x1fffff)+((r>>42)&0x1fffff))*_wynorm-3.0f;}

static	inline	void	wymum(uint64_t	*A,	uint64_t	*B){	uint64_t	hh=(*A>>32)*(*B>>32),	hl=(*A>>32)*(uint32_t)*B,	lh=(uint32_t)*A*(*B>>32),	ll=(uint64_t)(uint32_t)*A*(uint32_t)*B;	*A=((hl>>32)|(hl<<32))^hh;	*B=((lh>>32)|(lh<<32))^ll;}

static	inline	uint64_t	wyhash64(uint64_t	A,	uint64_t	B){	A^=0xa0761d6478bd642full;	B^=0xa0761d6478bd642full;	wymum(&A,&B);	A^=0xa0761d6478bd642full;	B^=0xa0761d6478bd642full;	wymum(&A,&B);	return	A^B;}

template<uint64_t	N>
static	inline	void	tiger(float	*w,	float	*g,	float	*m,	float	lr){
	float	r=lr*2;
	#pragma omp parallel for
	for(uint64_t	i=0;	i<N;	i++){	m[i]+=0.03125f*(g[i]-m[i]);	w[i]-=r*((m[i]>0)-0.5f);	}
}

template<uint64_t	N>
struct	Data{
	float	*data;
#ifndef		use_aligned_malloc
	Data(){	data=(float*)aligned_alloc(64,N*sizeof(float));	}
#else	
	Data(){	data=(float*)_aligned_malloc(N*sizeof(float),64);	}
#endif	
	~Data(){	free(data);	}
	void	save(FILE	*F){	fwrite(data,N*sizeof(float),1,F);	}
	bool	load(FILE	*F){	return	fread(data,N*sizeof(float),1,F)==1;	}
	void	zero(void){	memset(data,0,N*sizeof(float));	}
	void	rand(float	norm=1){	for(uint64_t	i=0;	i<N;	i++)	data[i]=norm*wy2gau(wyrand(&global.prng));	}
};

template<uint64_t	I,	uint64_t	O,	uint64_t	C>
struct	linear{
	static	Data<I*O>	g;
	Data<I*O>	w,m;
	linear(){	w.rand(1.0f/sqrtf(I));	m.zero();	}
	void	fw(Data<I*C>	&inp,	Data<O*C>	&out){
		Eigen::Map<Eigen::Matrix<float,I,O>,Eigen::Aligned64>	mw(w.data,I,O);
		Eigen::Map<Eigen::Matrix<float,I,C>,Eigen::Aligned64>	mi(inp.data,I,C);
		Eigen::Map<Eigen::Matrix<float,O,C>,Eigen::Aligned64>	mo(out.data,O,C);
		mo.noalias()=(1/sqrtf(I))*(mw.transpose()*mi);
	}
	
	void	bk(Data<I*C>	&inp,	Data<O*C>	&gin,	Data<I*C>	&gra,	bool	accumulate=false){
		Eigen::Map<Eigen::Matrix<float,I,O>,Eigen::Aligned64>	mw(w.data,I,O),	mg(g.data,I,O);
		Eigen::Map<Eigen::Matrix<float,I,C>,Eigen::Aligned64>	mi(inp.data,I,C),	gr(gra.data,I,C);
		Eigen::Map<Eigen::Matrix<float,O,C>,Eigen::Aligned64>	gi(gin.data,O,C);
		if(accumulate)	gr.noalias()+=(1/sqrtf(I))*(mw*gi);
		else	gr.noalias()=(1/sqrtf(I))*(mw*gi);
		mg.noalias()=(1/sqrtf(I*C))*(mi*gi.transpose());
		tiger<I*O>(w.data,g.data,m.data,global.learning_rate);
	}
};
template<uint64_t	I,	uint64_t	O,	uint64_t	C>
Data<I*O>	linear<I,O,C>::g;

}

#endif
