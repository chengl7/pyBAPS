import numpy as np
from distlink.common.server import Server
from multiprocessing import Process, Pipe
import socket
from distlink.worker.block import BlockProcess
from distlink.common.misc import disp_usage_forever

class LocalServer(Server):
    """LocalServer class, one per machine."""
    def __init__(self, parentConn=None, parentAddress=None, authKey=None, 
                 serverName=None, nChild=0, childBlockList=None, logFile=None):
        """Call server superclass initialization, start memMonitor, start BlockProcess children."
            
        Args:
            parentConn: connection to parent.
            parentAddress: address of parent.
            authKey: authentication key.
            serverAddress: address of this server.
            serverName: name of this server.
            nChild: number of child servers.
            childBlockList: blocks belonging to children.
            logFile: file to write logs.
        """

        super().__init__(parentConn=parentConn, parentAddress=parentAddress, authKey=authKey, 
             serverName=serverName, nChild=0, childBlockList=None, logFile=logFile)
        # get connections from children
        self.nChild = nChild
        self.childUpdateFlagArr = np.ones(nChild, dtype=bool)
        self.childConns = list()
        self.childBlockList = childBlockList   # a list indicating the blocks for each child
        
        if logFile:
            self.memMonitor = Process(target=disp_usage_forever,args=(self.logger.info,),name=socket.gethostname())
            self.memMonitor.start()
        else:
            self.memMonitor = None
        
        for i in range(nChild):
            parConn,chiConn=Pipe()
            self.childConns.append(parConn)
            blockIndex= childBlockList[i][0]
            serverName = f'block-process-({blockIndex[0]},{blockIndex[1]})'
            bp = BlockProcess(parentConn=chiConn, serverName=serverName, blockIndex=blockIndex)
            bp.start()
            
    def close(self):
        """Terminate memMonitor and close server."""
        if self.memMonitor:
            self.memMonitor.terminate()
        super().close()

