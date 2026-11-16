package krigingsTests;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.Map;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.hortonmachine.gears.io.shapefile.OmsShapefileFeatureReader;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.junit.Test;

import krigingsPointCase.Krigings;

/** New option doLocalVariogram of krigingsPointCase.Krigings: finite results, kriging variance, independence of the thread count. */
public class TestLocalVariogram {

	private Map<Integer, double[]> run(int threads, HashMap<Integer, double[]>[] var) throws Exception {
		OmsShapefileFeatureReader sr = new OmsShapefileFeatureReader();
		sr.file = "resources/Input/krigings/PointCase/rainstations.shp"; sr.readFeatureCollection();
		OmsShapefileFeatureReader pr = new OmsShapefileFeatureReader();
		pr.file = "resources/Input/krigings/PointCase/basins_passirio_width0.shp"; pr.readFeatureCollection();
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/krigings/PointCase/rain_test.csv"; reader.idfield = "ID"; reader.tStart = "2000-01-01 00:00";
		reader.tTimestep = 60; reader.fileNovalue = "-9999"; reader.initProcess();
		Krigings k = new Krigings();
		k.inStations = sr.geodata; k.fStationsid = "ID_PUNTI_M"; k.inInterpolate = pr.geodata; k.fInterpolateid = "netnum";
		k.doLocalVariogram = true; k.nVario = 10; k.nKrig = 5; k.threads = threads; k.boundedToZero = true;
		Map<Integer, double[]> all = new HashMap<>();
		int step = 0;
		while (reader.doProcess) {
			reader.nextRecord(); k.inData = reader.outData; k.executeKriging();
			for (Map.Entry<Integer, double[]> e : k.outData.entrySet()) all.put(step * 100000 + e.getKey(), e.getValue());
			if (var != null && step == 0) var[0] = k.outKrigingVariance;
			step++;
		}
		reader.close();
		return all;
	}

	@Test
	public void testLocalVariogramDeterministic() throws Exception {
		@SuppressWarnings("unchecked") HashMap<Integer, double[]>[] var = new HashMap[1];
		Map<Integer, double[]> a = run(1, var), b = run(4, null);
		assertTrue(a.size() > 0);
		for (Map.Entry<Integer, double[]> e : a.entrySet()) {
			double v = e.getValue()[0];
			assertTrue("finite estimate", Double.isFinite(v) && v >= 0.0);
			assertEquals(v, b.get(e.getKey())[0], 0.0);
		}
		for (double[] v : var[0].values()) assertTrue("kriging variance >= 0", v[0] >= 0.0);
	}
}
