//
///*!
// * \file      comms.c
// *
// * \brief     Comms subsytem functions
// *
// *
// *
// * \code
// *
// * 				 _______    ______    ____    ____    ____    ____     ______
// * 				/ ______)  /  __  \  |    \  /    |  |    \  /    |   / _____)
// * 			   / /         | |  | |  |  \  \/  /  |  |  \  \/  /  |  ( (____
// *            ( (          | |  | |  |  |\    /|  |  |  |\    /|  |   \____ \
// *             \ \______   | |__| |  |  | \__/ |  |  |  | \__/ |  |   _____) )
// *              \_______)  \______/  |__|      |__|  |__|      |__|  (______/
// *
// *
// * \endcode
// *
// * \author    Daniel Herencia
// *
// * \author    Robert Molina
// */
//
//#include <comms.h>
//
///*------IMPORTANT----------*/
///*
// * STATE CONTINGENCY ONLY RX
// * CHECK SLEEP MODE OF SX1262 WHILE NOT TX OR RX
// * MULTY THREAD
// * MATRIX OF READ-SALOMON => FADING IN SEVERAL PACKETS
// * SF AND CRC STORED IN MEMORY AND CHANGED BY TELECOMMANDS
// * send sf and crc at telemetry packet
// * THINK TELEMETRY PACKET!!!!!!!!!!
// * If stop sending in the middle of a window => set to zero count_packet[]
// * CAD_TIMER_TIMEOUT => CHANGE THE VALUE IN COMMS.H (and the next 6 definitions too)
// */
//
////#include "sx126x-hal.h"
//
//static RadioEvents_t RadioEvents;	//SHOULD THIS BE IN MAIN??? IS TO HANDLE IRQ???
//
//uint32_t air_time;
//uint8_t Buffer[BUFFER_SIZE];
//
///*
// * To avoid the variables being erased if a reset occurs, we have to store them in the Flash memory
// * Therefore, they have to be declared as a single-element array
// */
//uint8_t calib_packets = 0; //Counter of the calibration packets received
//uint8_t tle_packets = 0; //Counter of the tle packets received
//uint8_t telemetry_packets = 0; //Counter of telemetry packets sent
//
//uint8_t count_packet[] = {0};	//To count how many packets have been sent (maximum WINDOW_SIZE)
//uint8_t count_window[] = {0};	//To count the window number
//uint8_t count_rtx[] = {0};		//To count the number of retransmitted packets
//uint8_t i = 0;					//Auxiliar variable for loop
//uint8_t j=0;
//uint8_t k=0;
//
//uint64_t ack;					//Information rx in the ACK (FER DESPLAÇAMENTS DSBM)
//uint8_t nack_number;			//Number of the current packet to retransmit
//bool nack;						//True when retransmition necessary
//bool full_window;				//Stop & wait => to know when we reach the limit packet of the window
//bool statemach = true; //If true, comms workflow follows the statemachine. This value should be controlled by OBC
//							//Put true before activating the statemachine thread. Put false before ending comms thread
//bool send_data = false; //If true, the state machine send packets every airtime
//bool send_telemetry = false; //If true, we have to send telemetry packets instead of payload data
//uint8_t num_telemetry = 0; //Total of telemetry packets that have to be sent (computed when telecomand send Telemetry received)
//bool contingency = false; //True if we are in contingency state => only receive
//uint8_t uno=1; //TRUE 0x01
////uint8_t false=0;
//uint8_t SF=7;
////uint8_t Buffer[BUFFER_SIZE];
////bool PacketReceived = false;
////bool RxTimeoutTimerIrqFlag = false;
//
//
///*
// * STATE MACHINE VARIABLES
// */
//typedef enum
//{
//    LOWPOWER,      // low power consumption operating mode when we end always the transceiver is not transmitting nor processing received packets.
//    RX,      // CAD (channel activity detected) interruption occurs, call configuration()-->stateMachine(), START_CAD
//    RX_TIMEOUT, //transmitting process, trying to receive a packet a timer is activated so as to avoid getting suck in a bug or other receiving error.
//    RX_ERROR, //when an error in the receiving process occurs, we ruturn to START_CAD state.
//    TX,  //transmit data and telemetry to the ground station (tx_function() --> LOWPOWER), TxDone IRQ, when the sate is changed from the process_telecomand().
//    TX_TIMEOUT, //Tx Timmer produces a timeout exception. Happen when the transmission gets stuck or any error during the transmission occurs -->LOWPOWER
//    START_CAD, // enter in START_CAD state from any of the other reception states
//}States_t;
//
//typedef enum
//{
//    CAD_FAIL,
//    CAD_SUCCESS,
//    PENDING,
//}CadRx_t;
//
//States_t State = LOWPOWER;   //Variable to store the current state
//
//int8_t RssiValue = 0;		//RSSI computed value
//int8_t SnrValue = 0;		//SnrValue computed value
//
//CadRx_t CadRx = CAD_FAIL; 		//Current CAD state
//bool PacketReceived = false;		//To know if a packet have been received or not
//bool RxTimeoutTimerIrqFlag = false;	////Flag to know if there has been any interruption
//
//int16_t RssiMoy = 0;			//RSSI stored value
//int8_t SnrMoy = 0;			//SnrMoy stored value
//uint16_t RxCorrectCnt = 0;  ////Counter of correct received packets
//uint16_t BufferSize = BUFFER_SIZE; ////Buffer size
//
////TIMERS
//TimerEvent_t CADTimeoutTimer;
//TimerEvent_t RxAppTimeoutTimer;
//
///*
// *  ---------------- END ----------------
// */
//
//void configuration(void){
//
//	Radio.Init( &RadioEvents );
//
//	Radio.SetChannel( RF_FREQUENCY );
//
//	uint8_t sf;
//	Read_Flash(SF_ADDR, &sf, 1);
//	uint8_t coding_rate;
//	Read_Flash(CRC_ADDR, &coding_rate, 1);
//
//	Radio.SetTxConfig( MODEM_LORA, TX_OUTPUT_POWER, 0, LORA_BANDWIDTH,
//								   LORA_SPREADING_FACTOR, LORA_CODINGRATE,
//								   LORA_PREAMBLE_LENGTH, LORA_FIX_LENGTH_PAYLOAD_ON,
//								   true, 0, 0, LORA_IQ_INVERSION_ON, TX_TIMEOUT_VALUE );	//In the original example it was 3000
//
//
//	//SHALL WE CARE ABOUT THE RX TIMEOUT VALUE??? IF YES, CHANGE IT IN SetRx FUNCTION
//	Radio.SetRxConfig( MODEM_LORA, LORA_BANDWIDTH, LORA_SPREADING_FACTOR,
//								   LORA_CODINGRATE, 0, LORA_PREAMBLE_LENGTH,
//								   LORA_SYMBOL_TIMEOUT, LORA_FIX_LENGTH_PAYLOAD_ON,
//								   0, true, 0, 0, LORA_IQ_INVERSION_ON, true );
//
//
//	//Air time calculus
//	air_time = Radio.TimeOnAir( MODEM_LORA , PACKET_LENGTH );
//
//	Flash_Read_Data( COUNT_PACKET_ADDR, count_packet , sizeof(count_packet) );
//	//Read from Flash count_packet
//	Flash_Read_Data( COUNT_WINDOW_ADDR, count_window , sizeof(count_window) );
//	//Read from Flash count_window
//	Flash_Read_Data( COUNT_RTX_ADDR, count_rtx , sizeof(count_rtx) );
//	//Read from Flash count_rtx
//	ack = 0xFFFFFFFFFFFFFFFF; //Initially confifured 111..111
//	nack = false;
//	State = RX;
//
//
//};
//
////CAD: CHANNEL ACTIVITY DETECTED
//
//void tx_function(void){  //this function verify if there are packets left to be sent in this transmitting windows, and the calls the packaging function. It also updates the variables in the flash memory.
//	//configuration();
//	if (!full_window)
//	{
//		packaging(); //Start the TX by packaging all the data that will be transmitted
//		//SX126xSetPayload(); //Aquesta fa el writebuffer, sha de posar direccions com a la pag 48 del datasheet
//		Radio.Send( Buffer, BUFFER_SIZE );
//		Flash_Write_Data( COUNT_PACKET_ADDR , count_packet , sizeof(count_packet) ); //Read from Flash count_packet
//		Flash_Write_Data( COUNT_WINDOW_ADDR , count_window , sizeof(count_window) ); //Read from Flash count_window
//		Flash_Write_Data( COUNT_RTX_ADDR , count_rtx , sizeof(count_rtx) ); //Read from Flash count_rtx
//	}
//};
//
//
///* I THINK THAT THIS FUNCTION IS NOT NEEDED*/
//void rx_function(void){     //receives the packet detected
//	Radio.Rx( RX_TIMEOUT_VALUE );
//
//};
//
//void packaging(void){
//	//NACK packets at the beginnig of the next window
//
//	if (nack)
//	{
//		while(i<sizeof(ack))
//		{
//			//function to obtain the packets to retx
//			if(!((ack >> i) & 1)) //When position of the ack & 1 != 1 --> its a 0 --> NACK
//			{
//				nack_number = i;	//Current packet to rtx
//				//Packet from last window => count_window - 1
//				Flash_Read_Data( PHOTO_ADDR + (count_window[0]-1)*WINDOW_SIZE*BUFFER_SIZE + (nack_number)*BUFFER_SIZE , Buffer , sizeof(Buffer) );	//Direction in HEX
//				count_rtx[0]++;
//			}
//			i++;
//		}
//		if (i==sizeof(ack)){
//			i=0;
//			ack=0xFFFFFFFFFFFFFFFF;
//			nack = false;
//		}
//
//	}
//	else if (send_telemetry){
//		Flash_Read_Data( TELEMETRY_ADDR + telemetry_packets*(UPLINK_BUFFER_SIZE-1) , Buffer , sizeof(Buffer) );
//		Radio.Send( Buffer, BUFFER_SIZE );
//		telemetry_packets++;
//		if (telemetry_packets == num_telemetry){
//			send_telemetry = false;
//		}
//	}
//	else //no NACKS nor telemetry
//	{
//		Flash_Read_Data( PHOTO_ADDR + count_window[0]*WINDOW_SIZE*BUFFER_SIZE + (count_packet[0]-count_rtx[0])*BUFFER_SIZE , Buffer , sizeof(Buffer) );	//Direction in HEX
//		if (count_packet[0] < WINDOW_SIZE - 1)
//		{
//			count_packet[0]++;
//		}
//		else
//		{
//			count_packet[0] = 0;
//			count_window[0]++;
//			full_window = true;
//		}
//	}
//};
//
//
///*This function is called when a new photo is stored in the last photo position*/
//void resetCommsParams(void){
//	count_packet[0] = 0;
//	count_window[0] = 0;
//	count_rtx[0] 	= 0;
//}
//
//
//void stateMachine(void){
//    uint16_t PacketCnt = 0;
//
//    RadioEvents.TxDone = OnTxDone;
//    RadioEvents.RxDone = OnRxDone;
//    RadioEvents.TxTimeout = OnTxTimeout;
//    RadioEvents.RxTimeout = OnRxTimeout;
//    RadioEvents.RxError = OnRxError;
//    RadioEvents.CadDone = OnCadDone;
//
//    //Timer used to restart the CAD
//    TimerInit( &CADTimeoutTimer, CADTimeoutTimeoutIrq );
//    TimerSetValue( &CADTimeoutTimer, CAD_TIMER_TIMEOUT );
//
//    //App timmer used to check the RX's end
//    TimerInit( &RxAppTimeoutTimer, RxTimeoutTimerIrq );
//    TimerSetValue( &RxAppTimeoutTimer, RX_TIMER_TIMEOUT );
//
//    configuration();
//
//    SX126xConfigureCad( CAD_SYMBOL_NUM, CAD_DET_PEAK,CAD_DET_MIN, CAD_TIMEOUT_MS);      // Configure the CAD
//    Radio.StartCad( );                                                                  // do the config and lunch first CAD
//
//    State=RX;
//    switch(State)
//    {
//        case RX_TIMEOUT:
//        {
//			#if(FULL_DBG)
//            	printf( "RX Timeout\r\n");
//			#endif
//            //RxTimeoutCnt++;
//            State = START_CAD;
//         break;
//        }
//        case RX_ERROR:
//        {
//			#if(FULL_DBG)
//            	printf( "RX Error\r\n");
//			#endif
//            //RxErrorCnt++;
//            PacketReceived = false;
//            State = START_CAD;
//        break;
//        }
//        case RX:
//        {
//            if( PacketReceived == true )
//            {
//                PacketReceived = false;     // Reset flag
//                RxCorrectCnt++;         // Update RX counter
//                process_telecommand(Buffer[0], Buffer[1]); //We send the buffer to be used only in the cases of info = 1 byte
//				#if(FULL_DBG)
//					printf( "Rx Packet n %d\r\n", PacketCnt );
//				#endif
//                State = START_CAD;
//            }
//            else
//            {
//                if (CadRx == CAD_SUCCESS)
//                {
//                    //channelActivityDetectedCnt++;   // Update counter
//					#if(FULL_DBG)
//                    	printf( "Rxing\r\n");
//					#endif
//                    RxTimeoutTimerIrqFlag = false;
//                    TimerReset(&RxAppTimeoutTimer);	// Start the Rx's's Timer
//                    Radio.Rx( RX_TIMEOUT_VALUE );
//                }
//                else
//                {
//                    TimerStart(&CADTimeoutTimer);   // Start the CAD's Timer
//                }
//                State = LOWPOWER;
//            }
//            break;
//        }
//// mirar que passa si no arriba el packet(ACK)
//
//        case TX:
//        {
//            printf("Send Packet n %d \r\n",PacketCnt);
//            if( PacketCnt == 0xFFFF)
//            {
//                PacketCnt = 0;
//            }
//            else
//            {
//                PacketCnt ++;
//            }
//            //Send Frame
//            DelayMs( 1 );
//            if (send_data){
//            	tx_function();
//            }
//            State = LOWPOWER;
//            break;
//        }
//        case TX_TIMEOUT:
//        {
//            State = LOWPOWER;
//            break;
//        }
//        case START_CAD:
//        {
//            //i++;    // Update NbTryCnt
//            TimerStop(&RxAppTimeoutTimer);  // Stop the Rx's Timer
//            // Trace for debug
//            if(CadRx == CAD_FAIL)
//            {
//				#if(FULL_DBG)
//            		printf("No CAD detected\r\n");
//				#endif
//            }
//            CadRx = CAD_FAIL;           // Reset CAD flag
//            //DelayMs(randr(10,500));     //Add a random delay for the PER test => CHECK THIS WARNING
//            DelayMs(100);
//			#if(FULL_DBG)
//            	printf("CAD %d\r\n",i);
//			#endif
//            Radio.StartCad( );          //StartCad Again
//            State = LOWPOWER;
//        break;
//        }
//        case LOWPOWER:
//        default:
//        	//this is a low power consumption operating mode where we end always the transceiver is not transmitting nor processing received packets, and no IRQ are needed to be handled.
//        	//Then, the only way to exit this mode is changing the state manually from another function, or that an interruption request jumps from the transceiver.
//
//            break;
//    }
//}
//
//
//
//
///*
// * FUNCTIONS OBTAINED FROM EXAMPLE MAIN.C
// */
//
//
//void OnTxDone( void )
//{
//    Radio.Standby( );
//    State = TX;
//}
//
//void OnRxDone( uint8_t *payload, uint16_t size, int16_t rssi, int8_t snr )
//{
//    Radio.Standby( );
//    BufferSize = size;
//    memcpy( Buffer, payload, BufferSize );
//    RssiValue = rssi;
//    SnrValue = snr;
//    PacketReceived = true;
//    RssiMoy = (((RssiMoy * RxCorrectCnt) + RssiValue) / (RxCorrectCnt + 1));
//    SnrMoy = (((SnrMoy * RxCorrectCnt) + SnrValue) / (RxCorrectCnt + 1));
//    State = RX;
//}
//
//void OnTxTimeout( void )
//{
//    Radio.Standby( );
//    State = TX_TIMEOUT;
//}
//
//void OnRxTimeout( void )
//{
//    Radio.Standby( );
//    if( RxTimeoutTimerIrqFlag )
//    {
//        State = RX_TIMEOUT;
//    }
//    else
//    {
//		#if(FULL_DBG)
//        	printf(".");
//		#endif
//        Radio.Rx( RX_TIMEOUT_VALUE );   //  Restart Rx
//        //SymbTimeoutCnt++;               //  if we pass here because of Symbol Timeout
//        State = LOWPOWER;
//    }
//}
//
//void OnRxError( void )
//{
//    Radio.Standby( );
//    State = RX_ERROR;
//}
//
//void OnCadDone( bool channelActivityDetected)
//{
//    Radio.Standby( );
//
//    if( channelActivityDetected == true )
//    {
//        CadRx = CAD_SUCCESS;
//    }
//    else
//    {
//        CadRx = CAD_FAIL;
//    }
//    State = RX;
//}
//
//static void CADTimeoutTimeoutIrq( void )
//{
//    Radio.Standby( );
//    State = START_CAD;
//}
//
//static void RxTimeoutTimerIrq( void )
//{
//    RxTimeoutTimerIrqFlag = true;
//}
//
//void setContingency(bool cont){
//	contingency=cont;
//}
//void sendTelemetry(){
//	if(!contingency){
//		send_telemetry=false;
//		num_telemetry=(uint8_t) 34+BUFFER_SIZE + 1; //cast to integer to erase the decimal part
//		State = TX;
//	}
//	Write_Flash(TELEMETRY_ADDR, &uno , 1);
//}
//void sendData() {
//	if (!contingency){ //com entrar estat contingency
//	State = TX;
//	send_data = true;
//	}
//
//}
//void stopsendingData(){
//	send_data = false;
//	count_packet[0] = 0;
//
//}
//void sendcalibration(){
//	uint8_t calib[UPLINK_BUFFER_SIZE-1]; //RX
//	for (i=1; i<UPLINK_BUFFER_SIZE; i++){
//		calib[k-1]=Buffer[k];
//	}
//	Write_Flash(CALIBRATION_ADDR, calib, sizeof(calib));
//	calib_packets = calib_packets + 1;
//	uint8_t integer_part = (uint8_t) 138/UPLINK_BUFFER_SIZE;
//	if(calib_packets == integer_part+1){
//	calib_packets = 0;
//	}
//
//}
//void ackData(){
//	ack = ack & Buffer[1];
//	for(j=2; j<ACK_PAYLOAD_LENGTH; j++){
//	ack = (ack << 8*j) & Buffer[j];
//	}
//	count_window[0] = 0;
//	full_window = false;
//	if (ack != 0xFFFFFFFFFFFFFFFF){
//	nack = true;
//	}
//	State = TX;
//
//}
//int n=0;
//void setTime(){
//	uint8_t time[4];  //4bytes
//			for(n=0; n<4; n++){
//				time[n]=Buffer[n+1];
//			}
//			Write_Flash(SET_TIME_ADDR, time, sizeof(time));
//}
//void tle(){
//	/*Aquí si es reben els TLEs en diferents paquets doncs s'han d'anar
//			 * escrivint poc a poc o anar emmagatzemant la info i quan es tinguin
//			 * tots els bytes junts s'emmagatzemen (suposo que els paquets arriben
//			 * seguits)*/
//			uint8_t tle[UPLINK_BUFFER_SIZE-1];
//			for (k=1; k<UPLINK_BUFFER_SIZE; k++){
//			tle[k-1]=Buffer[k];
//			}
//			Write_Flash(TLE_ADDR + tle_packets*UPLINK_BUFFER_SIZE, tle, sizeof(tle));
//			tle_packets++;
//			uint8_t integer_part = (uint8_t) 138/UPLINK_BUFFER_SIZE;
//			if (tle_packets == integer_part+1){
//			tle_packets = 0;
//			}
//
//}
//void process_telecommand(uint8_t header, uint8_t info) {
//	switch(header) {
//	case RESET2:
//		/*Segons el drive s'ha de fer el reset si val 1 el bit, així que potser
//		 * s'hauria de posar un if*/
//		HAL_NVIC_SystemReset();
//		break;
//	case NOMINAL:
//		Write_Flash(NOMINAL_ADDR, &info, 1);
//		break;
//	case LOW:
//		Write_Flash(LOW_ADDR, &info, 1);
//		break;
//	case CRITICAL:
//		Write_Flash(CRITICAL_ADDR, &info, 1);
//		break;
//	case EXIT_LOW_POWER:
//		Write_Flash(EXIT_LOW_POWER_FLAG_ADDR, &info, 1);
//		Write_Flash(EXIT_LOW_ADDR, &uno, 1);
//		break;
//	case SET_TIME:
//		setTime();
//		break;
//	case SET_CONSTANT_KP:
//		Write_Flash(KP_ADDR, &info, 1);
//		break;
//	case TLE:
//		tle();
//		break;
//	case SET_GYRO_RES:
//		/*4 possibles estats,rebrem 00/01/10/11*/
//		Write_Flash(GYRO_RES_ADDR, &info, 1);
//		break;
//	case SENDDATA:
//		sendData();
//		break;
//	case SENDTELEMETRY:
//		sendTelemetry();
//		break;
//	case STOPSENDINGDATA:
//		stopsendingData();
//		break;
//	case ACKDATA:
//		ackData();
//		//translate the ACK/NACK into an understandable format for the transmission function, so that when the packaging function is called, it sends the packets from the last window that needs to be retransmitted. After changing the ACK format, switches the states to TX.
//		break;
//	case SET_SF:  //semicolon added in order to be able to declare SF here
//						//	/*4 cases (4/5, 4/6, 4/7,1/2), so we will receive and store 0, 1, 2 or 3*/
//		if (info == 0) SF = 7;
//		else if (info == 1) SF = 8;
//		else if (info == 2) SF = 9;
//		else if (info == 3) SF = 10;
//		else if (info == 4) SF = 11;
//		else if (info == 5) SF = 12;
//		Write_Flash(SF_ADDR, &SF, 1);
//		Write_Flash(CRC_ADDR, &Buffer[2], 1);
//		break;
//	case SEND_CALIBRATION:
//		sendcalibration();
//		break;
//	case TAKEPHOTO:
//		Write_Flash(PAYLOAD_STATE_ADDR, &uno, 1);
//		Write_Flash(PL_TIME_ADDR, &info, 4);
//		Write_Flash(PHOTO_RESOL_ADDR, &Buffer[5], 1);
//		Write_Flash(PHOTO_COMPRESSION_ADDR, &Buffer[6], 1);
//		break;
//	case TAKERF:
//		Write_Flash(PAYLOAD_STATE_ADDR, &uno, 1);
//		Write_Flash(PL_TIME_ADDR, &info, 8);
//		Write_Flash(PHOTO_RESOL_ADDR, &Buffer[5], 1);
//		Write_Flash(F_MIN_ADDR, &Buffer[9], 1);
//		Write_Flash(F_MAX_ADDR, &Buffer[10], 1);
//		Write_Flash(DELTA_F_ADDR, &Buffer[11], 1);
//		Write_Flash(INTEGRATION_TIME_ADDR, &Buffer[12], 1);
//		break;
//	case SEND_CONFIG: ; //semicolon added in order to be able to declare SF here
//		uint8_t config[CONFIG_SIZE];
//		Read_Flash(CONFIG_ADDR, config, CONFIG_SIZE);
//		Radio.Send(config, CONFIG_SIZE);
//		break;
//	}
//	//ha arribat telecomand (a obc), adressa o flag
//}
