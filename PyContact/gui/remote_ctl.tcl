#------------------------------------------------------------------
# $Id: remote_ctl.tcl,v 1.6 2003/02/12 21:33:11 oliver Exp $
# based on bounce.tcl and vmdcollab.tcl
# from http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/vmdcollab/
#
# start this in VMD and send commands to the listening port to have VMD
# execute them remotely
# (also see http://www.tcl.tk/scripting/netserver.html)
#
# Usage: vmd -e remote_ctl.tcl
# or vmd> source remote_ctl.tcl
#
# Security: we only allow connections from localhost (see acpt)
#
# Bugs:
# * once a wrong command was sent, the connection appears
# to 'block' and does not accept correct commands later
# * does not write result back to socket (one way connection...) so
# there is no way to inquire objects in vmd

namespace eval remote_ctl {
  variable main
  variable clients
  variable default_vmd_port
  set default_vmd_port 5555
  # I am too dumb to set the default value for port from
  # $default_vmd_port so I put 5555 in there literally
  proc start { {port 5050} } {
    variable main
    set main [socket -server remote_ctl::acpt $port]
    putlog "Listening on port $port"
  }
  proc acpt { sock addr port } {
    variable clients
    if {[string compare $addr "127.0.0.1"] != 0} {
      putlog "Unauthorized connection attempt from $addr port $port"
      close $sock
      return
    }
    putlog "Accept $sock from $addr port $port"
    set clients($sock) 1
    fconfigure $sock -buffering line
    fileevent $sock readable [list remote_ctl::recv $sock]
  }
  proc recv { sock } {
    variable main
    variable clients
    if { [eof $sock] || [catch {gets $sock line}]} {
      # end of file or abnormal connection drop:
      # shut down this connection
      close $sock
      putlog "Closing $sock"
      unset clients($sock)
      } else {
        if {[string compare $line "quit"] == 0} {
          # prevent new connections
          # existing connections stay open
          # No -- Bug(?): 'quit' closes VMD...
          putlog "Disallowing incoming connections by request of $sock"
          close $main
        }
        # execute the received commands
        # should check for runtime errors which otherwise leave the connection
        # in an unusable state
        #eval $line
        set rc [catch $line result]
        if { $rc } {
          #puts $sock "Error executing comand '$line': n$result"
          puts "Error executing comand '$line': $result"
          } else {
            #puts $sock $result
            puts $result
            if {[string first "set" $line] == 0} {

            }

          }
        }
      }
      ###### would like the last line from stdout in line ###########
      # (or any working solution....)
      proc send { sock line} {
        variable clients
        # send reply to connecting client
        putlog "send '$line' to $sock"
        puts $sock $line
      }

      proc putlog { text } {
        puts $text
        return
      }

      proc addToList { l numb } {
        set result {}
        foreach i $numb {
          lappend result [expr $i + $numb]
        }
        return $result
      }
    }

    remote_ctl::putlog "Starting remote_ctl server in vmd: connect with something like"
    remote_ctl::putlog "telnet localhost 5050"
    remote_ctl::start
