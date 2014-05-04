#include "adjsteptest.hh"

int ra_a_rand(RARR *arr)
{
    int len, dim, i;
    ra_ndim(arr, &dim);
    if ( dim <= 0 ) return 0;

    len = 1;
    for ( i=0; i<dim; ++i ) len *= arr->_dims[i].n0;
    for ( i=0; i<len; ++i ) {
        arr->_s0[i] = 1.0*rand()/RAND_MAX * 10.0;
    }
    return 0;
}

int ra_a_copy(const RARR *arr, RARR *arr1)
{
    int len, dim, dim1, i;
    ra_ndim(arr,  &dim);
    ra_ndim(arr1, &dim1);

    if (dim == dim1){
        len=1;
        for ( i=0; i<dim; ++i ) len *= arr->_dims[i].n0;
        for ( i=0; i<len; ++i ) arr1->_s0[i] =  arr->_s0[i];
        return 0;
    }
    return 1;
}

int ra_a_swap(RARR *arr, RARR *arr1){
    int len, dim, dim1, i;
    ireal tmp;
    ra_ndim(arr,  &dim);
    ra_ndim(arr1, &dim1);
    
    if (dim == dim1){
        len=1;
        for ( i=0; i<dim; ++i ) len *= arr->_dims[i].n0;
        for ( i=0; i<len; ++i ) {
            tmp = arr1->_s0[i];
            arr1->_s0[i] =  arr->_s0[i];
            arr->_s0[i] = tmp;
        }
        return 0;
    }
    return 1;
    
}

void iwave_rdom_rand(IWAVE *state)
{
    int i;
    for (i=0; i<(state->model).ld_a.narr; i++){
        ra_a_rand(&((state->model).ld_a._s[i]));
    }
}

void rdom_rand(RDOM *rd)
{
    int i;
    for (i=0; i<rd->narr; i++){
        ra_a_rand(&(rd->_s[i]));
    }
}

void iwave_rdom_copy(const IWAVE *state, IWAVE *state1)
{
    int i;
    for (i=0; i<(state->model).ld_a.narr; i++)
        if(ra_a_copy(&((state->model).ld_a._s[i]), &((state1->model).ld_a._s[i])))
            cout << "ERROR: two RARR have different dim!! cannot copy.\n";
}

void rdom_copy(const RDOM *rd_src, RDOM *rd_dst)
{
    int i;
    for (i=0; i<(rd_src->narr); i++)
        if(ra_a_copy(&(rd_src->_s[i]), &(rd_dst->_s[i])))
            cout << "ERROR: two RARR have different dim!! cannot copy.\n";
}

void adjrelation1(RDOM * iwf, RDOM * iwa,
                 std::vector<IWAVE> &iwtmp, ostream & str){
    
    RARR *upd_n = &(iwf->_s[D_UC]);
    RARR *ucd_n = &(iwf->_s[D_UP]);
    RARR *upb_o = &(iwtmp[1].model.ld_c._s[D_UP]);
    RARR *ucb_o = &(iwtmp[1].model.ld_c._s[D_UC]);
    
    ireal ip, ip1;
    ra_a_inner(upd_n, upb_o, &ip);
    ra_a_inner(ucd_n, ucb_o, &ip1);
    
    ireal axnorm, ynorm, tmp;
    ra_a_inner(upd_n, upd_n, &axnorm);
    ra_a_inner(ucd_n, ucd_n, &tmp);
    axnorm = sqrt(axnorm+tmp);
    ra_a_inner(upb_o, upb_o, &ynorm);
    ra_a_inner(ucb_o, ucb_o, &tmp);
    ynorm = sqrt(ynorm+tmp);
    
    RARR *upd_o  = &(iwtmp[0].model.ld_c._s[D_UP]);
    RARR *ucd_o  = &(iwtmp[0].model.ld_c._s[D_UC]);
    RARR *csqd_o = &(iwtmp[0].model.ld_c._s[D_CSQ]);
    RARR *ucb_n  = &(iwa->_s[D_UC]);
    RARR *upb_n  = &(iwa->_s[D_UP]);
    RARR *csqb_n = &(iwa->_s[D_CSQ]);
    
    ireal ip2, ip3, ip4;
    ra_a_inner(upd_o, upb_n, &ip2);
    ra_a_inner(csqd_o, csqb_n, &ip3);
    ra_a_inner(ucd_o, ucb_n, &ip4);
    
    str << "< Ax,    y > = " << ip + ip1<< endl;
    str << "<  x, A^Ty > = " << ip2 + ip3 + ip4<< endl;
    str << " |Ax|  = " << axnorm << endl;
    str << "  |y|  = " << ynorm << endl;
    str << "< Ax, y> - < x, A^Ty> \n";
    str << "--------------------- = " << abs(-ip -ip1 + ip2 + ip3 + ip4)/axnorm/ynorm << endl;
    str << "       |Ax||y|        \n";
    str << "100*machine eps = " << 100 * numeric_limits<float>::epsilon() << endl;
}

void adjrelation2(std::vector<RDOM *> iwf,
                  std::vector<RDOM *> iwa,
                  std::vector<IWAVE> &iwtmp, ostream & str){
    
    if ((iwf.size()+iwa.size()) != iwtmp.size() ) {
        str << "input IWAVE vectors have wrong dimensions:\n";
        str << "iwf.size()   = " << iwf.size() << endl;
        str << "iwa.size()   = " << iwa.size() << endl;
        str << "iwtmp.size() = " << iwtmp.size() << endl;
        str << "RIGHT relation: iwf.size()+iwa.size()=iwtmp.size()!\n";
    }
        
    // Ax
    RARR *ucdd_n = &(iwf[1]->_s[D_UP]);
    RARR *updd_n = &(iwf[1]->_s[D_UC]);
    // y
    RARR *ucdb_o = &(iwtmp[3].model.ld_c._s[D_UC]);
    RARR *updb_o = &(iwtmp[3].model.ld_c._s[D_UP]);
    
    // < Ax,y >
    ireal ip, tmp1, tmp2;
    ra_a_inner(ucdd_n, ucdb_o, &ip);
    ra_a_inner(updd_n, updb_o, &tmp1);
    ip = ip + tmp1;
    
    // norm
    ireal axnorm, ynorm;
    ra_a_inner(ucdd_n, ucdd_n, &axnorm);
    ra_a_inner(updd_n, updd_n, &tmp1);
    axnorm = sqrt(axnorm+tmp1);
    ra_a_inner(ucdb_o, ucdb_o, &ynorm);
    ra_a_inner(updb_o, updb_o, &tmp1);
    ynorm = sqrt(ynorm+tmp1);
    
    // x
    RARR *ucdd_o  = &(iwtmp[2].model.ld_c._s[D_UC]);
    RARR *updd_o  = &(iwtmp[2].model.ld_c._s[D_UP]);
    RARR *csqdd_o = &(iwtmp[2].model.ld_c._s[D_CSQ]);
    RARR *ucd0_o  = &(iwtmp[0].model.ld_c._s[D_UC]);
    RARR *csqd0_o = &(iwtmp[0].model.ld_c._s[D_CSQ]);
        
    // A^Ty
    RARR *ucdb_n  = &(iwa[1]->_s[D_UC]);
    RARR *updb_n  = &(iwa[1]->_s[D_UP]);
    RARR *csqdb_n = &(iwa[1]->_s[D_CSQ]);
    RARR *ucb_n   = &(iwa[0]->_s[D_UC]);
    RARR *csqb_n  = &(iwa[0]->_s[D_CSQ]);
    
    // < x, A^Ty >
    ireal ip2, ip3;
    ra_a_inner(ucdd_o, ucdb_n, &ip2);
    ra_a_inner(updd_o, updb_n, &tmp1);
    ra_a_inner(csqdd_o, csqdb_n, &tmp2);
    ip2 = ip2  + tmp2 + tmp1;
    ra_a_inner(ucd0_o, ucb_n, &ip3);
    ra_a_inner(csqd0_o, csqb_n, &tmp2);
    ip3 = ip3 + tmp2;

    str << "< Ax,    y > = " << ip << endl;
    str << "<  x, A^Ty > = " << ip2 + ip3<< endl;
    str << " |Ax|  = " << axnorm << endl;
    str << "  |y|  = " << ynorm << endl;
    str << "< Ax, y> - < x, A^Ty> \n";
    str << "--------------------- = " << abs(-ip + ip2 +ip3)/axnorm/ynorm << endl;
    str << "       |Ax||y|        \n";
    str << "100*machine eps = " << 100 * numeric_limits<float>::epsilon() << endl;
}

void adjsteptest(std::vector<RDOM *> &iwf, std::vector<RDOM *> &iwa, IWaveInfo const &ic,
		  void (*tsf)(std::vector<RDOM *> , bool, int, void *fdpars), int iv,  
                  void *fdpars, PARARRAY * pars, FILE * stream){
    // int err;
    int n=iwf.size();
    // decleare and construct working space
    std::vector<IWAVE> iwtmp(n);
    for (int i=0; i<n; i++) {
        iwave_construct(&iwtmp[i],pars,stream,ic);
    }
    if(n == 2){
        // initialize iwa[1], D_CSQ to zero for adjoint operator
        ra_a_zero(&(iwa[n-1]->_s[D_CSQ]));
        rdom_copy(iwf[n-1], &(iwtmp[0].model.ld_c));
        rdom_copy(iwa[n-1], &(iwtmp[1].model.ld_c));
        
        if(ra_a_swap(&(iwa[n-1]->_s[D_UC]),&(iwa[n-1]->_s[D_UP]))!=0) cout << " ERROR!\n";
    }
    else if(n == 4){
        // initialize iwa[1] & iwa[3], D_CSQ to zero for adjoint operator
        ra_a_zero(&(iwa[n-1]->_s[D_CSQ]));
        ra_a_zero(&(iwa[1]->_s[D_CSQ]));
        ra_a_zero(&(iwa[1]->_s[D_UC]));
        rdom_copy(iwf[1], &(iwtmp[0].model.ld_c));
        rdom_copy(iwa[1], &(iwtmp[1].model.ld_c));
        rdom_copy(iwf[n-1], &(iwtmp[2].model.ld_c));
        rdom_copy(iwa[n-1], &(iwtmp[3].model.ld_c));
        
        if(ra_a_swap(&(iwa[1]->_s[D_UC]),&(iwa[1]->_s[D_UP]))!=0) cout << " ERROR in adjsteptest!\n";
        if(ra_a_swap(&(iwa[n-1]->_s[D_UC]),&(iwa[n-1]->_s[D_UP]))!=0) cout << " ERROR in adjsteptest!\n";
    }
    
    // call forward and adjoint operator
    tsf(iwf, true,  iv,  fdpars);
    tsf(iwa, false, iv,  fdpars);
    if (n>1){
    cout <<"====================================\n";
    cout <<"=========== adjoint test ===========\n";
    cout <<"====================================\n";
    if (n==2) {
        adjrelation1(iwf[n-1], iwa[n-1], iwtmp, cout);
    }
    else if (n==4) {
        std::vector<RDOM *> iwf0(2);
        std::vector<RDOM *> iwa0(2);
        iwf0.at(0) = iwf.at(1); iwf0.at(1) = iwf.at(3);
        iwa0.at(0) = iwa.at(1); iwa0.at(1) = iwa.at(3);
        adjrelation2(iwf0, iwa0, iwtmp, cout);
    }
    cout <<"========= end adjoint test =========\n";
    }

    for (int i=0; i<n; i++) {
      iwave_destroy(&iwtmp[i],ic.get_mdest());
    }
}


