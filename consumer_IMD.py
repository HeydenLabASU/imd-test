#!/usr/bin/env python3

import struct
import socket
from enum import Enum
import numpy as np
import time
import MDAnalysis as mda

u=mda.Universe("GMX/run.tpr","GMX/struct.gro")
stride=10

#All constants that are not originally defined in GROMACS source code
IMDV3 = 3
TRICLINIC_DIMENSIONS = 6
ORTHORHOMBIC_DIMENSIONS = 3

#Define constants
HEADERSIZE = 8
IMDVERSION = 2
ENERGYSIZE = 40

#Enum for IMDType
class IMDType(Enum):
    IMD_DISCONNECT = 0
    IMD_ENERGIES = 1
    IMD_FCOORDS = 2
    IMD_GO = 3
    IMD_HANDSHAKE = 4
    IMD_KILL = 5
    IMD_MDCOMM = 6
    IMD_PAUSE= 7
    IMD_TRATE = 8
    IMD_IOERROR = 9
    IMD_SESSIONINFO = 10
    IMD_RESUME = 11
    IMD_TIME = 12
    IMD_BOX = 13
    IMD_VELOCITIES = 14
    IMD_FORCES = 15

class IMDheader:
    def __init__(self, header_type, length):
        self.type = header_type
        self.length = length

    def pack(self):
        return struct.pack('!ii', self.type, self.length)

    def unpack(data):
        '''read header data from bytearray (network to host conversion)'''
        header_type, length = struct.unpack('!ii', data)
        return IMDheader(header_type, length)

    def unpack_le(data,header):
        '''read header.length from handshake (=IMDVERSION) as little endian'''
        length = struct.unpack_from('<i', data,4)
        return IMDheader(header.type, length[0])

    def unpack_be(data,header):
        '''read header.length from handshake (=IMDVERSION) as big endian'''
        length = struct.unpack_from('>i', data,4)
        return IMDheader(header.type, length[0])

class IMDEnergies:
    """
    typedef struct {
        int32 tstep;      /**< integer timestep index                    */
        float T;          /**< Temperature in degrees Kelvin             */
        float Etot;       /**< Total energy, in Kcal/mol                 */
        float Epot;       /**< Potential energy, in Kcal/mol             */
        float Evdw;       /**< Van der Waals energy, in Kcal/mol         */
        float Eelec;      /**< Electrostatic energy, in Kcal/mol         */
        float Ebond;      /**< Bond energy, Kcal/mol                     */
        float Eangle;     /**< Angle energy, Kcal/mol                    */
        float Edihe;      /**< Dihedral energy, Kcal/mol                 */
        float Eimpr;      /**< Improper energy, Kcal/mol                 */
    } IMDEnergies;
    """
    def __init__(self, tstep, T, Etot, Epot, Evdw, Eelec, Ebond, Eangle, Edihe, Eimpr):
        self.tstep = tstep
        self.T = T
        self.Etot = Etot
        self.Epot = Epot
        self.Evdw = Evdw
        self.Eelec = Eelec
        self.Ebond = Ebond
        self.Eangle = Eangle
        self.Edihe = Edihe
        self.Eimpr = Eimpr

    def unpack_be(data):
        '''read energy data in big endian format from bytearray'''
        variables = struct.unpack('>ifffffffff', data)
        return IMDEnergies(*variables)

    def unpack_le(data):
        '''read energy data in little endian format from bytearray'''
        variables = struct.unpack('<ifffffffff', data)
        return IMDEnergies(*variables)

class IMDfcoords:
    def __init__(self, fcoords, natoms):
        self.fcoords = fcoords
        self.natoms = natoms

    def unpack_le(u,data,n):
        '''read fcoords in little endian format from bytearray'''
        #print(f'receiving {n} coordinates/forces')
        #fcoords = struct.iter_unpack("<fff",data)
        #for a in u.atoms:
        #    a.position = struct.unpack("<fff",data)
        format="<"+n*"fff"
        #fcoords = np.asarray(struct.unpack(format,data), dtype=float)
        #u.atoms.positions = fcoords.reshape(n,3)
        u.trajectory.ts._pos[:] = np.reshape(struct.unpack(format,data),(n,3)) #line profiler
        fcoords = u.trajectory.ts._pos
        return IMDfcoords(fcoords, n)

    def unpack_be(u,data,n):
        '''read fcoords in big endian format from bytearray'''
        #print(f'receiving {n} coordinates/forces')
        #fcoords = struct.iter_unpack(">fff",data)
        #for a in u.atoms:
        #    a.position = struct.unpack(">fff",data)
        format=">"+n*"fff"
        #fcoords = np.asarray(struct.unpack(format,data), dtype=float)
        #u.atoms.positions = fcoords.reshape(n,3)
        u.trajectory.ts._pos[:] = np.reshape(struct.unpack(format,data),(n,3))
        fcoords = u.trajectory.ts._pos
        return IMDfcoords(fcoords, n)

class IMDBox:
    def __init__(self, box_dimensions):
        self.box_dimensions = box_dimensions

    def unpack_le(data, num_elements):
        '''read box dimensions in little endian format from bytearray'''
        box_dimensions = struct.unpack(f'<{num_elements}f', data)
        return IMDBox(box_dimensions)

    def unpack_be(data, num_elements):
        '''read box dimensions in big endian format from bytearray'''
        box_dimensions = struct.unpack(f'>{num_elements}f', data)
        return IMDBox(box_dimensions)

def imd_readn(sock, n):
    """Read n bytes from a socket."""
    buffer = bytearray(n)
    view = memoryview(buffer)
    nleft = n

    while nleft > 0:
        chunk = sock.recv(nleft)
        if not chunk:
            break
        view[:len(chunk)] = chunk
        view = view[len(chunk):]
        nleft -= len(chunk)

    return bytes(buffer)

def imd_writen(sock, data):
    """Write n bytes to a socket."""
    nleft = len(data)
    total_sent = 0

    while nleft > 0:
        sent = sock.send(data[total_sent:])
        if sent == 0:
            raise RuntimeError("Socket connection broken")
        total_sent += sent
        nleft -= sent

    return total_sent
   
def imd_send_pause(sock):
    """Send a pause message."""
    header = IMDheader(IMDType.IMD_PAUSE.value, 0)
    packed_header = header.pack()
    sock.sendall(packed_header)
    print("IMD_PAUSE signal sent")
    #fill_header(header, IMDType.IMD_PAUSE.value, 0)
    #return imd_writen(sock, header.pack()) == HEADERSIZE

def imd_send_go(sock):
    header = IMDheader(IMDType.IMD_GO.value, 0)
    packed_header = header.pack()
    sock.sendall(packed_header)
    print("IMD_GO signal sent")

def imd_send_trate(sock,stride):
    header = IMDheader(IMDType.IMD_TRATE.value, stride)
    packed_header = header.pack()
    sock.sendall(packed_header)
    print(f'IMD_TRATE = {stride} sent')

def handle_client_connection(sock,stride):
    print("Client connected")
    paused = False

    #Receive and interpret handshake
    header_data = imd_readn(sock, HEADERSIZE)
    if not header_data:
        print("Failed to receive handshake header")
        return

    endian=-1
    header = IMDheader.unpack(header_data)
    if header.type == IMDType.IMD_HANDSHAKE.value:
        print(f'Received handshake')
        header = IMDheader.unpack_le(header_data,header)
        if(header.length == IMDV3):
            endian=0
            print(f'little endian detected')
        else:
            header = IMDheader.unpack_be(header_data,header)
            if(header.length == IMDV3):
                endian=1
                print(f'big endian detected')
            else:
                print(f'ERROR: could not detect endianess')
        
        if(endian!=-1):
            print("sending IMD_GO signal")
            imd_send_go(sock)
    
    if(stride!=0):
        imd_send_trate(sock,stride)

    with mda.Writer("imd-test.trr", len(u.atoms)) as w:    
        while True:
            #Receive header
            header_data = imd_readn(sock, HEADERSIZE)
            if not header_data:
                print("Failed to receive header")
                break
            header = IMDheader.unpack(header_data)
            if header.type == IMDType.IMD_DISCONNECT.value:
                print("Client requested disconnect")
                break
            elif header.type == IMDType.IMD_ENERGIES.value:
                energy_data = imd_readn(sock, ENERGYSIZE)
                if(endian == 0):
                    energy = IMDEnergies.unpack_le(energy_data)
                else:
                    energy = IMDEnergies.unpack_be(energy_data)
                #print("Energy: ", energy.tstep, energy.T, energy.Etot, energy.Epot)
            elif header.type == IMDType.IMD_FCOORDS.value:
                fcoords_data = imd_readn(sock, header.length*12)
                if(endian == 0):
                    fcoords = IMDfcoords.unpack_le(u,fcoords_data, header.length)
                else:
                    fcoords = IMDfcoords.unpack_be(u,fcoords_data, header.length)
                #print(f'{header.length}\n')
                #for x in fcoords.fcoords:
                #    print(f'X {x[0]} {x[1]} {x[2]}')
                w.write(u)
            elif header.type == IMDType.IMD_BOX.value:
                print(f'{header.length}')
                box_dimensions_data = imd_readn(sock, header.length * 4)
                if header.length == TRICLINIC_DIMENSIONS:
                    print('triclinic box')
                elif header.length == ORTHORHOMBIC_DIMENSIONS:
                    print('orthorhombic box')

                if endian == 0:
                    box_dimensions = IMDBox.unpack_le(box_dimensions_data, header.length)
                else:
                    box_dimensions = IMDBox.unpack_be(box_dimensions_data, header.length)        
                print(f'{box_dimensions.box_dimensions}')
            else:
                print("IMD message with unknown header.type received")
                print(f'message type = {header.type}')
                break
    sock.close()
    print("Disconnected") 

def run_consumer():
    consumer_address = ('localhost', 8888)
    consumer_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    consumer_socket.connect(consumer_address)
    time.sleep(1)
    handle_client_connection(consumer_socket,stride)

if __name__ == "__main__":
    run_consumer()
