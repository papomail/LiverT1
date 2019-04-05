classdef patients < handle
  
   properties (SetAccess = private)
      HN
      ID
      
   end
  
   events
      InsufficientFunds
   end
   methods
      function PD = LoadPatient(accNum,initBal)
         PD.HN = accNum;
         PD.ID = initBal;
         PD.Niftis = 
        % BA.AccountListener =  AccountManager.addAccount(BA);
      end
      function deposit(PD,amt)
         PD.AccountBalance = PD.AccountBalance + amt;
         if PD.AccountBalance > 0
            PD.AccountStatus = 'open';
         end
      end
      function withdraw(PD,amt)
         if (strcmp(PD.AccountStatus,'closed')&& PD.AccountBalance <= 0)
            disp(['Account ',num2str(PD.AccountNumber),' has been closed.'])
            return
         end
         newbal = PD.AccountBalance - amt;
         PD.AccountBalance = newbal;
         if newbal < 0
            notify(PD,'InsufficientFunds')
         end
      end
      function getStatement(PD)
         disp('-------------------------')
         disp(['Account: ',num2str(PD.AccountNumber)])
         ab = sprintf('%0.2f',PD.AccountBalance);
         disp(['CurrentBalance: ',ab])
         disp(['Account Status: ',PD.AccountStatus])
         disp('-------------------------')
      end
   end
   
end
