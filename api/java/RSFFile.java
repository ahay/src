package rsf;

public abstract class RSFFile{

    public static final int MAX_DIMS = m8rConstants.SF_MAX_DIM;
    
    public SWIGTYPE_p_sf_File file; //do not access unless you know what you are doing
    
    public RSFFile(String filename, boolean input){
    /* create a new rsf file */
        if (input) {
            file = m8r.sf_input(filename);
        } else {
            file = m8r.sf_output(filename);
        }
    }
    
    public void close(){
            m8r.sf_fileclose(this.file);
    }
         
    public int getN(int index){
            String key = "n"+Integer.toString(index);
            return getInt(key);
    }

    public void setLabel(int index, String label){
            String key = "label" + Integer.toString(index);
            setString(key,label);
    }
    
    public String getLabel(int index){
            String key = "label" + Integer.toString(index);
            return getString(key);
    }

    public void setUnit(int index, String unit){
            String key = "unit"+Integer.toString(index);
            setString(key,unit);
    }
    
    public String getUnit(int index){
            String key = "unit" + Integer.toString(index);
            return getString(key);
    }
    
    public void setDelta(int index, float delta){
            String key = "d"+Integer.toString(index);
            setFloat(key,delta);
    }
    
    public float getDelta(int index){
            String key = "d"+Integer.toString(index);
            return getFloat(key);
    }

    public void setOrigin(int index, float origin){
            String key = "o"+Integer.toString(index);
            setFloat(key,origin);
    }
    
    public float getOrigin(int index){
            String key = "o"+Integer.toString(index);
            return getFloat(key);
    }

    public void setN(int index, int n){
            String key = "n"+Integer.toString(index);
            setInt(key,n);
    }
    
    private void setInt(String key, int n){
            try {
                m8r.sf_putint(this.file,key,n);
            } catch (Exception e){
                System.err.println(e);
            }
    }
    
    private void setFloat(String key, float n){
        try {
            m8r.sf_putfloat(this.file,key,n);
        } catch (Exception e){
            System.err.println(e);
        }
    }
    
    private void setString(String key, String val){
        try {
            m8r.sf_putstring(this.file,key,val);
        } catch (Exception e){
            System.err.println(e);
        }
    }
    
    private int getInt(String key){
        try {
            int[] value = new int[]{0};
            boolean found = m8r.sf_histint(this.file,key,value);
            if (found) return value[0];
            else return 1;
        } catch (Exception e){
            System.err.println(e);
            return 0;
        }
    }
    
    private float getFloat(String key){
        try {
            float[] value = new float[]{0.0f};
            boolean found = m8r.sf_histfloat(this.file,key,value);
            if (found) return value[0];
            else return 0.0f;
        } catch (Exception e){
            System.err.println(e);
            return 0.0f;
        }
    }
    
    private String getString(String key){
        try {
            String value = m8r.sf_histstring(this.file,key);
            if (value == null) return "";
            else return value;
        } catch (Exception e){
            System.err.println(e);
            return null;
        }   
    }
}
