package net.sf.osra;

import javax.management.RuntimeErrorException;

/**
 * Pure JNI bridge for OSRA library.
 * 
 * @author <a href="mailto:dmitry.katsubo@gmail.com">Dmitry Katsubo</a>
 */
class OsraLibJni extends OsraLib {

	private static final String	NAME	= "osra_java";

	static {
		try {
			System.loadLibrary(NAME);
		}
		catch (UnsatisfiedLinkError e) {
			throw new RuntimeErrorException(e, "Check that lib" + NAME + ".so/" + NAME
						+ ".dll is in PATH or in java.library.path (" + System.getProperty("java.library.path") + ")");
		}
	}
}
