select 
  run.triggerconfig,
  run.brtimestamp as date,
  run.runnumber,
  run.runtype,
  run.updateunixtime-run.brunixtime as time_in_run,
  run.eventsinrun,
  trigger.name,
  trigger.enabled,
  trigger.scaledown
from 
  run,trigger 
where 
  run.runnumber=trigger.runnumber 
  and trigger.name = 'CLOCK' 
  and run.runtype = 'VERNIERSCAN' 
  and run.eventsinrun > 7000000 
  and (run.triggerconfig like '%Run12%' or run.triggerconfig like '%Run15%') 
order by 
  run.runnumber desc limit 100;

