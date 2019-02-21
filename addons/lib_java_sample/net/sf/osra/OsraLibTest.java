package net.sf.osra;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;

import org.apache.commons.io.IOUtils;
import org.junit.Test;

/**
 * Sample usage of the library.
 * 
 * @author <a href="mailto:dmitry.katsubo@gmail.com">Dmitry Katsubo</a>
 */
public class OsraLibTest {

	@Test
	public void testProcessImage() throws IOException {
		StringWriter writer = new StringWriter();
		InputStream is = new BufferedInputStream(new FileInputStream("test/test.png"));

		byte[] imageData = IOUtils.toByteArray(is);

		int result = OsraLibJni.processImage(imageData, writer, 0, false, 0, 0, 0, false, false, "sdf", "inchi", true,
					true, true, true, true);

		System.out.println("OSRA completed with result:" + result + " structure:\n" + writer.toString() + "\n");
	}
}
