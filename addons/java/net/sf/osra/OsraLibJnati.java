package net.sf.osra;

import java.io.IOException;
import java.io.InputStream;
import java.util.PropertyResourceBundle;

import net.sf.jnati.NativeCodeException;
import net.sf.jnati.deploy.NativeLibraryLoader;

/**
 * JNI bridge for OSRA library based on JNATI library.
 * 
 * @author <a href="mailto:dmitry.katsubo@gmail.com">Dmitry Katsubo</a>
 */
public class OsraLibJnati extends OsraLibJni {

	private static final String	NAME	= "osra";

	private static final String	VERSION;

	public static String getVersion() {
		return VERSION;
	}

	static {
		try {
			VERSION = getVersionFromResource();

			NativeLibraryLoader.loadLibrary(NAME, VERSION);
		}
		catch (NativeCodeException e) {
			// Unable to handle this.
			throw new RuntimeException(e);
		}
		catch (IOException e) {
			// Unable to handle this.
			throw new RuntimeException(e);
		}
	}

	private static final String	MAVEN_PROPERTIES	= "META-INF/maven/net.sf.osra/osra/pom.properties";

	/**
	 * This function looks up the OSRA version from Maven properties file.
	 */
	private static String getVersionFromResource() throws IOException {
		InputStream is = OsraLibJnati.class.getClassLoader().getResourceAsStream(MAVEN_PROPERTIES);

		try {
			PropertyResourceBundle resourceBundle = new PropertyResourceBundle(is);

			return resourceBundle.getString("version");
		}
		finally {
			is.close();
		}
	}
}
