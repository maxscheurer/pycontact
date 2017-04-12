/*CPP Class for chTimer.h*/

#include <chTimer.h>

class ChTimer {

public:
	ChTimer(){};
	~ChTimer(){};
	
	//
	// Start the Timer
	//
	int 
	start() { chTimerGetTime( &m_start ); return 0; };
	
	//
	// Stop the Timer
	//
	int 
	stop()  { chTimerGetTime( &m_end ); return 0; };
	
	//
	// Get elapsed Time
	//
	double 
	getTime() { return chTimerElapsedTime( &m_start, &m_end ); };
	double 
	getTime(int iterations) { return (chTimerElapsedTime( &m_start, &m_end )/static_cast<double>(iterations)); };
	//
	// Get Bandwidth
	//
	double 
	getBandwidth(double size) { return chTimerBandwidth( &m_start, &m_end, size); };
	double 
	getBandwidth(double size, int iterations) { return (chTimerBandwidth( &m_start, &m_end, size)*static_cast<double>(iterations)); };
	

private:
	chTimerTimestamp m_start;
	chTimerTimestamp m_end;

};
