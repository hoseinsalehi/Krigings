package pokriging.engine;

import java.lang.reflect.Field;

/**
 * Command-line runner without OMS (tests and batch use): java pokriging.engine.Cli key=value ...
 * Keys are the fields of {@link KrigingConfig} (= the @In keys of pokriging.PoKriging).
 */
public final class Cli {
    public static void main(String[] args) throws Exception {
        KrigingConfig c = new KrigingConfig();
        for (String a : args) {
            int k = a.indexOf('=');
            set(c, a.substring(0, k), a.substring(k + 1));
        }
        KrigingEngine e = new KrigingEngine(c);
        e.run();
        System.out.println("done: " + e.info);
    }

    static void set(Object o, String key, String v) throws Exception {
        Field f = o.getClass().getField(key);
        Class<?> t = f.getType();
        if (t == int.class) f.setInt(o, Integer.parseInt(v));
        else if (t == double.class) f.setDouble(o, Double.parseDouble(v));
        else if (t == boolean.class) f.setBoolean(o, Boolean.parseBoolean(v));
        else f.set(o, v);
    }
}
