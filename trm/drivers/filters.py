#!/usr/bin/env python

"""
Routines to talk to the ULTRASPEC filter wheel

-- really need to have wheel to interact with it.

"""

import pyserial

class FilterWheel (object):

    # usual port names
    PORT_NAMES = [ 
        '/dev/filterwheel', # Linux
        '/dev/tty.usbserial-FTTJI800', # OS X
        '/dev/ttyUSB0', # another name
        ]

    # Seconds to block while waiting for port open
    TIME_OUT = 2

    # Default bits per second for COM port.
    DATA_RATE = 19200

    def __init__(self):
        self.connected  = False
        self.response   = ''
        self.ID         = ''
        self.serialPort = None

    def connect(self):
        # default connection: 
        self.serialPort = serial.Serial(0, DATA_RATE, timeout=TIME_OUT)
			throw new WheelException("Error: Could not find Serial port.");
		}	

		try {
			// open serial port, and use class name for the appName.
			serialPort = (SerialPort) portId.open(this.getClass().getName(),
					TIME_OUT);

			if(serialPort == null) throw new WheelException("port not opened");
			// set port parameters
			serialPort.setSerialPortParams(DATA_RATE,
					SerialPort.DATABITS_8,
					SerialPort.STOPBITS_1,
					SerialPort.PARITY_NONE);

			// open the streams
			input = serialPort.getInputStream();
			//input = new BufferedReader(new InputStreamReader(serialPort.getInputStream()));
			output = serialPort.getOutputStream();
					
				
		} catch (PortInUseException e){
		    System.err.println("Error: port in use");
		    throw new WheelException("Error: port in use");
		} catch (Exception e) {
		    System.err.println(e.toString());
		    throw e;
		}
		connected=true;
	}

    SerialPort serialPort

	/** Buffered input stream from the port */
	private InputStream input;
	/** The output stream to the port */
	private OutputStream output;





    public synchronized boolean isConnected(){
	return connected;
    }

	/**
	 * This should be called when you stop using the port.
	 * This will prevent port locking on platforms like Linux.
	 */
	public synchronized void close() throws Exception{
		try{
			this.sendCommand("WEXITS");
		}catch(Exception e){
			System.out.println("warning: may have failed to stop serial mode on filter wheel");
		}
		if (serialPort != null) {
			try{
				serialPort.close();
				input.close();
				output.close();
				connected=false;
			}catch (Exception e){
				System.err.println("Error: cannot close port");
				throw new WheelException("Error: cannot close port");
			}
		}
	}

	// Sends command to the wheel, and awaits a response for up to 30s- this is blocking
	public synchronized String sendCommand(String comm) throws Exception{
		String res = "";
		if (serialPort != null) {
			// first strip all trailing whitespace off command
			comm = comm.replaceAll("/n","");
			
			//determine appropriate delay (ms)
			int delay = 30000;
			if(comm.equals("HOME") || comm.equals("GOTO")){
				delay = 30000;
			}
			
			// add CR/LF
			comm += "\0x0D\0x0A";
			try{
				// send string
				output.write(comm.getBytes());
				output.flush();
							
				// need this for some reason, don't know why
				Thread.sleep(100);
			
				// wait up to 1s for a response
				int i=0;
				while(i<100){
					int numB = input.available();
					if(numB > 0) {
						break;
					}
					i++;
					// sleep for delay/100
					Thread.sleep(delay/100);
				}
				if(i==100){
				    // strip characters off end
				    comm = comm.replaceAll("\0x0D\0x0A","");
				    throw new WheelException("Command '"+comm+"' timed out");
				}
			
				int numBytes = input.available();
				byte[] buffer = new byte[numBytes];
				input.read(buffer,0,numBytes);
				res = new String(buffer);
				// strip CR/LF from response
				res = res.replaceAll("\n","");
				res = res.replaceAll("\r","");
			
			}catch(Exception e){
				System.err.println(e.toString());
				throw new WheelException(e.toString());
			}
		}else{
			throw new WheelException("Filter wheel not connected via serial");
		}
		return res;
	}
	
	//enables Wheel for serial responses
	public synchronized void init() throws Exception{
		response = this.sendCommand("WSMODE");
		if(!response.equals("!")){
			throw new WheelException("Couldn't initialise wheel for serial commands");
		}
	}
	
	//homes filter wheel, returns ID of wheel
	public synchronized void home() throws Exception{
		response = this.sendCommand("WHOME");
		if(response.equals("ER=3")){
			throw new WheelException("Could not ID filter wheel after HOME");
		}	
		if(response.equals("ER=1")){
			throw new WheelException("Filter wheel home took too many steps");
		}	
		ID = response;
	}
	
	// returns letter ID of current filter wheel
	public synchronized String getID() throws Exception{
		response = this.sendCommand("WIDENT");
		if (!response.matches("A|B|C|D|E")){
			throw new WheelException("Bad filter wheel ID: "+response);
		}
		return response;
	}
	
	// gets wheel position (from 1 to 6)
	public synchronized int getPos() throws Exception{
		response = this.sendCommand("WFILTR");
		int filtNum = new Integer(response);
		return filtNum;
	}
	
	// get filter wheel names - these should always be 1-6
	public synchronized int[] getNames() throws Exception{
		response = this.sendCommand("WREAD");
		String[] tmp = response.split("\\s+");
		int[] retArr = new int[tmp.length];
		for(int i=0; i< tmp.length; i++){
			retArr[i] = new Integer(tmp[i]);
		}
		return retArr;
	}
	
	// move filter wheel to desired location
	public synchronized void move(int pos) throws Exception{
		response = this.sendCommand("WGOTO"+pos);
		if(response.equals("ER=5")){
			throw new WheelException("requested filter position ("+pos+") not valid");
		}
		if(response.equals("ER=4")){
			throw new WheelException("filter wheel is stuck");
		}
		if(response.equals("ER=6")){
			throw new WheelException("filter wheel is slipping");
		}
		if(!response.equals("*")){
			throw new WheelException();
		}
	}
	
	// filter wheel not responding? This should fix it
	public synchronized void reSpawn() throws Exception{
		try{
			this.close();
		}catch (Exception e) {
			System.out.println(e.toString());
		}
		Thread.sleep(1000);
		this.connect();
		this.init();
		this.home();
	}
	
	public static void main(String[] args) throws Exception {
		FilterWheel wheel = new FilterWheel();
		try{
			wheel.connect();
			wheel.init();
			System.out.println("Started");
			System.out.println(wheel.getPos());
			System.out.println(wheel.getID());
			int[] ids = wheel.getNames();
			for(int i=0; i<ids.length; i++){
				System.out.println("ID: "+ids[i]);
			}
			System.out.println("GOTO 9");
			try{
				wheel.move(9);
			}catch(Exception e){
				System.err.println(e.toString());
			}
			wheel.reSpawn();
			System.out.println("GOTO 2");
			try{
				wheel.move(2);
			}catch(Exception e){
				System.err.println(e.toString());
			}
			wheel.close();
			System.out.println("Finished normally");
			System.exit(0);						
		}catch( Exception e){
			System.out.println(e.toString());
			System.out.println("shutting down");
			wheel.close();
			System.exit(1);
		}
	}
}


class WheelException extends Exception
{
  String mistake;
//----------------------------------------------
// Default constructor - initializes instance variable to unknown
  public WheelException()
  {
    super();             // call superclass constructor
    mistake = "unknown";
  }
  
//-----------------------------------------------
// Constructor receives some kind of message that is saved in an instance variable.
  public WheelException(String err)
  {
   
    super(err);     // call super class constructor
    mistake = err;  // save message
  }
  
//------------------------------------------------  
// public method, callable by exception catcher. It returns the error message.
  public String getError()
  {
    return mistake;
  }
}


