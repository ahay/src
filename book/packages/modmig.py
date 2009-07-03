def modmig(par):
    Flow('layermodel',None,'vam nz=%{nz}d nx=%{nx}d dz=%{dz}g dx=%{dx}g' % par)
    Result('layermodel','grey 
