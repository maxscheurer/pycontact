/*
    Copyright (C) 2012  Aaron S. Keys

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* ***************************
	vmstream.h
	A lightweight C++ library to open a TCP connection to VMD 
	Aaron Keys, June 10, 2009

   ***************************/	

/** \mainpage VMD Streams Documentation

\section intro_sec Introduction 
The vmdstreams library provides a simple interface for opening a TCP connection 
to VMD from C++.  Once the TCP connection is established, the user can remote-
control VMD from a C++ code to automate analysis, control movies, or make custom
animations.  All that is required is some knowledge of VMDs Tcl interface (see 
http://www.ks.uiuc.edu/Research/vmd/current/ug/node104.html).

\section examples Examples

The vmdstream object works the same as an ostream in C++.  (That is, it works 
the same way as cout, cerr and ofstreams). Here is an example of drawing 20
spheres on the screen:

\code

#include <vmstream/vmdstream.h>

using namespace std;

int main( int argc, char *argv[] ) 
{ 
	//open a TCP connection to VMD
	vmdsock_t vmdsock = newvmdsock();
	vmdstream vmd(vmdsock); 
	
	//tell VMD to draw some spheres
	cout << "drawing spheres\n";
	
	for (int i=0; i<10; i++) {
		vmd << "draw sphere {" << i << " 0 0}" << endl;
	}

	for (int i=0; i<10; i++) {
		vmd << "draw sphere {0 " << i << " 0}" << endl;
	}

	//close the tcp connection
	cout << "finished" << endl; 
	closevmdsock(vmdsock);	
} 

\endcode

Notice that the "draw sphere" command is part of the VMD tcl scripting 
interface.  In fact, we can control almost every aspect of VMD using Tcl 
commands.

\section Installation

The vmdstreams codes consists of a single header file.  

\section References

This codes is based on the socket iostream library by David Ryan, released into 
public domain in Feb 2004.
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <sstream>
#include <fstream>
#include <string>
#include <streambuf> 
#include <iostream> 

#include <netdb.h> 
#include <sys/socket.h> 
#include <sys/ioctl.h> 
#include <netinet/in.h> 

#ifndef TIOCINQ 
#define TIOCINQ FIONREAD 
#endif 

class socket_buf : public std::streambuf 
{ 
private: 
   int _socket; 
   char _inBuffer[1024];  
   char _outBuffer[1024]; 

public: 
   socket_buf( int socket ) 
   { 
      _socket = socket; 
      setg( _inBuffer, _inBuffer, _inBuffer ); 
      setp( _outBuffer, _outBuffer+1024 ); 
   } 

protected: 
   int overflow(int x) 
   { 
#ifdef DEBUG_MODE
      std::cout << "socket_buf::overflow()" << x << std::endl; 
#endif
      return x; 
   } 
    
   int sync() 
   { 
#ifdef DEBUG_MODE
      std::cout << "socket_buf::sync()" << std::endl; 
#endif       
      // if no data availalbe just return. 
      if ( pbase() == pptr() ) 
         return 0; 
       
      // try and send the data. 
      int len = pptr() - pbase(); 
      int rc = send( _socket, pbase(), len, 0 ); 
      if ( rc < 0 ) 
         return rc; 
          
      setp( _outBuffer, _outBuffer+1024 ); 
      return 0; 
   } 

   int underflow() 
   { 
#ifdef DEBUG_MODE
      std::cout << "socket_buf::underflow() " << std::endl; 
#endif       
      if (gptr () < egptr ()) return *(unsigned char*)gptr (); 

      int len; 
       
      // find out how much data is available. 
      if( ioctl( _socket, TIOCINQ, &len) < 0 ) 
      { 
         // error 
         std::cerr << "socket_buf error" << std::endl; 
         len = 1; 
      } 
       
      // make sure length is atleast 1. We will block. 
      if ( len == 0 ) 
         len = 1; 
       
      // try and read in some data. 
      int read = recv( _socket, &_inBuffer, len, 0 ); 
      if ( read > 0 ) 
      { 
         setg( _inBuffer, _inBuffer, _inBuffer+read ); 
      } 
      else 
      { 
         return EOF;    
      } 
       
      return *(unsigned char*)gptr ();       
   } 

}; 


class siostream: public std::iostream 
{ 
private: 
   socket_buf _buf; 

public: 
   siostream( int socket ) 
   :_buf( socket ), std::iostream( &_buf ) 
   { 
   } 
}; 

#endif

typedef siostream vmdstream;

void writevmdinitscript(int port=5555, std::string extra_commands="")
{
	std::ofstream os("remote_ctl.tcl");
	
	os << "#------------------------------------------------------------------" << "\n";
	os << "# $Id: remote_ctl.tcl,v 1.6 2003/02/12 21:33:11 oliver Exp $" << "\n";
	os << "# based on bounce.tcl and vmdcollab.tcl" << "\n";
	os << "# from http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/vmdcollab/" << "\n";
	os << "#" << "\n";
	os << "# start this in VMD and send commands to the listening port to have VMD" << "\n";
	os << "# execute them remotely" << "\n";
	os << "# (also see http://www.tcl.tk/scripting/netserver.html)" << "\n";
	os << "#" << "\n";
	os << "# Usage: vmd -e remote_ctl.tcl" << "\n";
	os << "# or vmd> source remote_ctl.tcl" << "\n";
	os << "#" << "\n";
	os << "# Security: we only allow connections from localhost (see acpt)" << "\n";
	os << "#" << "\n";
	os << "# Bugs:" << "\n";
	os << "# * once a wrong command was sent, the connection appears" << "\n";
	os << "# to 'block' and does not accept correct commands later" << "\n";
	os << "# * does not write result back to socket (one way connection...) so" << "\n";
	os << "# there is no way to inquire objects in vmd" << "\n";
	os << "" << "\n";
	os << "namespace eval remote_ctl {" << "\n";
	os << "variable main" << "\n";
	os << "variable clients" << "\n";
	os << "variable default_vmd_port" << "\n";
	os << "set default_vmd_port 5555" << "\n";
	os << "# I am too dumb to set the default value for port from" << "\n";
	os << "# $default_vmd_port so I put 5555 in there literally" << "\n";
	os << "proc start { {port " << port << "} } {" << "\n";
	os << "variable main" << "\n";
	os << "set main [socket -server remote_ctl::acpt $port]" << "\n";
	os << "putlog \"Listening on port $port\"" << "\n";
	os << "}" << "\n";
	os << "proc acpt { sock addr port } {" << "\n";
	os << "variable clients" << "\n";
	os << "if {[string compare $addr \"127.0.0.1\"] != 0} {" << "\n";
	os << "putlog \"Unauthorized connection attempt from $addr port $port\"" << "\n";
	os << "close $sock" << "\n";
	os << "return" << "\n";
	os << "}" << "\n";
	os << "putlog \"Accept $sock from $addr port $port\"" << "\n";
	os << "set clients($sock) 1" << "\n";
	os << "fconfigure $sock -buffering line" << "\n";
	os << "fileevent $sock readable [list remote_ctl::recv $sock]" << "\n";
	os << "}" << "\n";
	os << "proc recv { sock } {" << "\n";
	os << "variable main" << "\n";
	os << "variable clients" << "\n";
	os << "if { [eof $sock] || [catch {gets $sock line}]} {" << "\n";
	os << "# end of file or abnormal connection drop:" << "\n";
	os << "# shut down this connection" << "\n";
	os << "close $sock" << "\n";
	os << "putlog \"Closing $sock\"" << "\n";
	os << "unset clients($sock)" << "\n";
	os << "} else {" << "\n";
	os << "if {[string compare $line \"quit\"] == 0} {" << "\n";
	os << "# prevent new connections" << "\n";
	os << "# existing connections stay open" << "\n";
	os << "# No -- Bug(?): 'quit' closes VMD..." << "\n";
	os << "putlog \"Disallowing incoming connections by request of $sock\"" << "\n";
	os << "close $main" << "\n";
	os << "}" << "\n";
	os << "# execute the received commands" << "\n";
	os << "# should check for runtime errors which otherwise leave the connection" << "\n";
	os << "# in an unusable state" << "\n";
	os << "# eval $line" << "\n";
	os << "set rc [catch $line result]" << "\n";
	os << "if { $rc } {" << "\n";
	os << "#puts $sock \"Error executing comand '$line': n$result\"" << "\n";
	os << "puts \"Error executing comand '$line': n$result\"" << "\n";
	os << "} else {" << "\n";
	os << "#puts $sock $result" << "\n";
	os << "#puts $result" << "\n";
	os << "}" << "\n";
	os << "}" << "\n";
	os << "}" << "\n";
	os << "###### would like the last line from stdout in line ###########" << "\n";
	os << "# (or any working solution....)" << "\n";
	os << "proc send { sock line} {" << "\n";
	os << "variable clients" << "\n";
	os << "# send reply to connecting client" << "\n";
	os << "putlog \"send '$line' to $sock\"" << "\n";
	os << "puts $sock $line" << "\n";
	os << "}" << "\n";
	os << "" << "\n";
	os << "proc putlog { text } {" << "\n";
	os << "puts $text" << "\n";
	os << "return" << "\n";
	os << "}" << "\n";
	os << "}" << "\n";
	os << "" << "\n";
	os << "remote_ctl::putlog \"Starting remote_ctl server in vmd: connect with something like\"" << "\n";
	os << "remote_ctl::putlog \"telnet localhost " << port << "\"" << "\n";
	os << "remote_ctl::start" << "\n";

	os << extra_commands;
	
	os.close();
}

typedef int vmdsock_t;

/**
\brief opens a TCP connection to VMD and returns the socket file descriptor
\param vmd_executable is the name of the vmd program as it is called from 
	the command line.  Typically, if vmd does not live in a global search path
	such as /usr/bin or /usr/local/bin, we make a dynamic link to the executable
	an place it in the /usr/bin directory.  For example, on mac os, the VMD 
	executable typically lives in an application package, for example:
	/Applications/VMD\ 1.8.6.app/Contents/vmd/vmd_MACOSXX86.  In this case,
	it is easiest to make a dynamic link to the executable using "ln" and place
	the link in the /usr/bin folder
\param port is the port to connect over
\param extra_commands is a string containing a list of Tcl commands to parse
	on startup
\return the socket file descriptor
*/
vmdsock_t newvmdsock(
	const char* vmd_executable="vmd", 
	int port=5555, 
	std::string extra_commands="")
{  
	writevmdinitscript(port, extra_commands);
	int rv, s;
        int count;

	pid_t  pid;
	pid = fork();
	if (pid == -1) {
		std::cerr << "ERROR: failed to fork child process\n";
		exit(1);
	}
	else if (pid == 0) {	//child	
		std::ostringstream command;
		command << vmd_executable <<" -e remote_ctl.tcl";
		int good = system(command.str().c_str());
		if (good != 0) {
			std::cerr << "ERROR: vmd executable failed to launch\n";
			std::cerr << "Please check your path: " << vmd_executable << "\n";
			exit(1);
		}
		std::cout << "VMD killed via quit command\n";
		exit(0);
	}	
	else {
		//sleep(2);	
		std::cout << "test socket" << std::endl; 

		struct hostent *he=gethostbyname("localhost"); 
		struct sockaddr_in sa; 
		memset(&sa, 0, sizeof(sa)); 
		sa.sin_family = PF_INET; 
		sa.sin_addr = *((struct in_addr *)he->h_addr);        
		sa.sin_port = htons(port); 
			 

                for (count=0, rv=1; rv && count<20; count++) {
		  s = socket(AF_INET, SOCK_STREAM, 0); 
		  rv = connect( s, (struct sockaddr*) &sa, sizeof(sa) ); 

		  std::cout << "connect = " << rv << std::endl;
		
		  if (rv) {
			  std::cerr << "TCP connection to vmd failed\n";
                          sleep(2);
                          close(s);
		  }
                  sleep(1);
                }
	}
	
	return s; 
}

/**
\brief closes the VMD socket connection
\param socketfd is the file descriptor of the socket to close
*/
void closevmdsock(vmdsock_t socketfd)
{
    shutdown(socketfd, SHUT_RDWR);
    close(socketfd);
}
