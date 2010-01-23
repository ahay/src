package rsf;

public class RSFHeader{

    public static final int MAX_DIMS = 3;
    
    public RSFHeader(){
        ndims = 0;
        esize = 4;
        path = "";
        format = "";
        deltas = new float[MAX_DIMS];
        origins = new float[MAX_DIMS];
        n = new int[MAX_DIMS];
        labels = new String[MAX_DIMS];
        units = new String[MAX_DIMS];
        name = "";
    }

    public void setName(String name){
            this.name = name;
    }
    

    public void setDelta(int index, float delta){
            try{
                    deltas[index-1] = delta;
            } catch (Exception e){
                    System.out.printf("Failed putting %f at index %d into deltas\n",delta,index);
                    System.exit(1);
            }
    }

    public void setOrigin(int index, float origin){
            try{
                    origins[index-1] = origin;
            } catch (Exception e){
                    System.out.printf("Failed putting %f at index %d into origins\n",origin,index);
            }
    }

    public void setN(int index, int n){
            try{
                    this.n[index-1] = n;
                    ndims += 1;
            } catch (Exception e){
                    System.out.printf("Failed putting %d at index %d into n\n",n,index);
            }
    }

    public void setLabel(int index, String label){
            try{
                    labels[index-1] = label;
            } catch (Exception e){
                    System.out.printf("Failed to put %s at index %d in labels\n",label,index);
            }
    }

    public void setUnit(int index, String unit){
            try{
                    units[index-1] = unit;
            } catch (Exception e){
                    System.out.printf("Failed to put %s at index %d in unit\n",unit,index);
            }
    }

    public void setPath(String path){
            this.path = path;
    }
    public void setFormat(String format){
            this.format = format;
    }

    public String toString(){
		String s = "dims: "+ndims+"\n"+intString("numbers",n)+floatString("origins",origins)+floatString("deltas",deltas)+ stringString("labels",labels)+stringString("units",units)+"path: " + path+"\n" + "format: "+ format+"\n"+"esize: "+esize+"\n";
		return s;
	}
	public String floatString(String name, float[] f){
		String s = name + ": [";
		for(int i = 0; i < f.length; ++i){
			s = s + String.valueOf(f[i])+" ";
		}
		s = s + "]\n";
		return s;
	}	
	public String intString(String name, int[] f){
		String s = name + ": [";
		for(int i = 0; i < f.length; ++i){
			s = s + String.valueOf(f[i])+" ";
		}
		s = s + "]\n";
		return s;
	}
	public String stringString(String name, String[] f){
		String s = name + ": [";
		for(int i = 0; i < f.length; ++i){
			s = s + f[i]+" ";
		}
		s = s + "]\n";
		return s;
	}

    public float[] deltas;
    public float[] origins;
    public int[] n;
    public String[] labels;
    public String[] units;
    public String path;
    public String format;
    public String name;
    public int esize;
    public int ndims;
}
