package rsf;

import java.util.Hashtable;

public class Par{

    public Par(String[] args){
        dictionary = new Hashtable<String,String>();
        for (String arg: args){
            try{
                if (arg.contains("=")){
                    String[] split = arg.split("=");
                    String key = split[0];
                    String val = split[1];
                    dictionary.put(key,val);
                    //System.out.printf("Added key: %s, value: %s\n",key,val);
                }
            } catch (Exception e){
                //do nothing
            }
        }
    }

    public boolean getBool(String key, boolean fallback){
        if (dictionary.containsKey(key)){
            Boolean t = new Boolean(fallback);
            try {
                String s = (String)dictionary.get(key);
                if (s.equals("y")){
                    return true;
                } else {
                    return false;
                }
            } catch (Exception e){
                //do nothing
            }
            return t.booleanValue();
        } else {
            return fallback;
        }
    }

    public int getInt(String key, int fallback){
        if(dictionary.containsKey(key)){
            Integer i = new Integer(fallback);
            try{
                String s = (String)dictionary.get(key);
                i = Integer.valueOf(s);
            } catch (Exception e){
                //do nothing
            }
            return i.intValue();
        } else {
            return fallback;
        }
    }

    public float getFloat(String key, float fallback){
        if(dictionary.containsKey(key)){
            Float f = new Float(fallback);
            try {
                String s = (String)dictionary.get(key);
                f = Float.valueOf(s);
            } catch (Exception e){
            }
            return f.floatValue();
        } else {
            return fallback;
        }
    }

    public String getString(String key, String fallback){
        if(dictionary.containsKey(key)){
           String s = fallback; 
           try {
                s = (String)dictionary.get(key);
           } catch (Exception e){
                //do nothing
           }
           return s;
        } else {
            return fallback;
        }
    }
    private Hashtable<String,String> dictionary;
}
