import torchtuples as tt
import torch
from torch import Tensor
import numpy as np
import pandas as pd

class SurvBase(tt.Model):
    """Base class for survival models. 
    Essentially same as torchtuples.Model, 
    """
    def predict_surv(self, input, batch_size=8224, numpy=None, eval_=True,
                     to_cpu=False, num_workers=0):
        """Predict the survival function for `input`.
        See `prediction_surv_df` to return a DataFrame instead.

        Arguments:
            input {dataloader, tuple, np.ndarray, or torch.tensor} -- Input to net.

        Keyword Arguments:
            batch_size {int} -- Batch size (default: {8224})
            numpy {bool} -- 'False' gives tensor, 'True' gives numpy, and None give same as input
                (default: {None})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            to_cpu {bool} -- For larger data sets we need to move the results to cpu
                (default: {False})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            [TupleTree, np.ndarray or tensor] -- Predictions
        """
        raise NotImplementedError

    def predict_surv_df(self, input, batch_size=8224, eval_=True, num_workers=0):
        """Predict the survival function for `input` and return as a pandas DataFrame.
        See `predict_surv` to return tensor or np.array instead.

        Arguments:
            input {dataloader, tuple, np.ndarray, or torch.tensor} -- Input to net.

        Keyword Arguments:
            batch_size {int} -- Batch size (default: {8224})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            pd.DataFrame -- Predictions
        """
        raise NotImplementedError

    def predict_hazard(self, input, batch_size=8224, numpy=None, eval_=True, to_cpu=False,
                       num_workers=0):
        """Predict the hazard function for `input`.

        Arguments:
            input {dataloader, tuple, np.ndarray, or torch.tensor} -- Input to net.

        Keyword Arguments:
            batch_size {int} -- Batch size (default: {8224})
            numpy {bool} -- 'False' gives tensor, 'True' gives numpy, and None give same as input
                (default: {None})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            grads {bool} -- If gradients should be computed (default: {False})
            to_cpu {bool} -- For larger data sets we need to move the results to cpu
                (default: {False})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            [np.ndarray or tensor] -- Predicted hazards
        """
        raise NotImplementedError

    def predict_pmf(self, input, batch_size=8224, numpy=None, eval_=True, to_cpu=False,
                    num_workers=0):
        """Predict the probability mass function (PMF) for `input`.

        Arguments:
            input {dataloader, tuple, np.ndarray, or torch.tensor} -- Input to net.

        Keyword Arguments:
            batch_size {int} -- Batch size (default: {8224})
            numpy {bool} -- 'False' gives tensor, 'True' gives numpy, and None give same as input
                (default: {None})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            grads {bool} -- If gradients should be computed (default: {False})
            to_cpu {bool} -- For larger data sets we need to move the results to cpu
                (default: {False})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            [np.ndarray or tensor] -- Predictions
        """
        raise NotImplementedError

class _CoxBase(SurvBase):
    
    duration_col = 'duration'
    event_col = 'event'
    truncation_col = 'truncation'

    def fit(self, input, target, batch_size=256, epochs=1, callbacks=None, verbose=True,
            num_workers=0, shuffle=True, metrics=None, val_data=None, val_batch_size=8224,
            **kwargs):
        """Fit  model with inputs and targets. Where 'input' is the covariates, and
        'target' is a tuple with (durations, events).
        
        Arguments:
            input {np.array, tensor or tuple} -- Input x passed to net.
            target {np.array, tensor or tuple} -- Target [durations, events]. 
        
        Keyword Arguments:
            batch_size {int} -- Elements in each batch (default: {256})
            epochs {int} -- Number of epochs (default: {1})
            callbacks {list} -- list of callbacks (default: {None})
            verbose {bool} -- Print progress (default: {True})
            num_workers {int} -- Number of workers used in the dataloader (default: {0})
            shuffle {bool} -- If we should shuffle the order of the dataset (default: {True})
            **kwargs are passed to 'make_dataloader' method.
    
        Returns:
            TrainingLogger -- Training log
        """
        self.training_data = tt.tuplefy(input, target)
        return super().fit(input, target, batch_size, epochs, callbacks, verbose,
                           num_workers, shuffle, metrics, val_data, val_batch_size,
                           **kwargs)

    def _compute_baseline_hazards(self, input, df, max_duration, batch_size, eval_=True, num_workers=0):
        raise NotImplementedError

    def target_to_df(self, target):
        durations, events, truncation = tt.tuplefy(target).to_numpy()
        df = pd.DataFrame({self.duration_col: durations, self.event_col: events, self.truncation_col: truncation}) 
        return df

    def compute_baseline_hazards(self, input=None, target=None, max_duration=None, sample=None, batch_size=8224,
                                set_hazards=True, eval_=True, num_workers=0):
        """Computes the Breslow estimates form the data defined by `input` and `target`
        (if `None` use training data).

        Typically call
        model.compute_baseline_hazards() after fitting.
        
        Keyword Arguments:
            input  -- Input data (train input) (default: {None})
            target  -- Target data (train target) (default: {None})
            max_duration {float} -- Don't compute estimates for duration higher (default: {None})
            sample {float or int} -- Compute estimates of subsample of data (default: {None})
            batch_size {int} -- Batch size (default: {8224})
            set_hazards {bool} -- Set hazards in model object, or just return hazards. (default: {True})
        
        Returns:
            pd.Series -- Pandas series with baseline hazards.
        """
        if (input is None) and (target is None):
            if not hasattr(self, 'training_data'):
                raise ValueError("Need to give a 'input' and 'target' to this function.")
            input, target = self.training_data
        df = self.target_to_df(target)#.sort_values(self.duration_col)
        if sample is not None:
            if sample >= 1:
                df = df.sample(n=sample)
            else:
                df = df.sample(frac=sample)
        input = tt.tuplefy(input).to_numpy().iloc[df.index.values]
        base_haz = self._compute_baseline_hazards(input, df, max_duration, batch_size,
                                                  eval_=eval_, num_workers=num_workers)
        if set_hazards:
            self.compute_baseline_cumulative_hazards(set_hazards=True, baseline_hazards_=base_haz)
        return base_haz

    def compute_baseline_cumulative_hazards(self, input=None, target=None, max_duration=None, sample=None,
                                            batch_size=8224, set_hazards=True, baseline_hazards_=None,
                                            eval_=True, num_workers=0):
        """See `compute_baseline_hazards. This is the cumulative version."""
        if ((input is not None) or (target is not None)) and (baseline_hazards_ is not None):
            raise ValueError("'input', 'target' and 'baseline_hazards_' can not both be different from 'None'.")
        if baseline_hazards_ is None:
            baseline_hazards_ = self.compute_baseline_hazards(input, target, max_duration, sample, batch_size,
                                                             set_hazards=False, eval_=eval_, num_workers=num_workers)
        assert baseline_hazards_.index.is_monotonic_increasing,\
            'Need index of baseline_hazards_ to be monotonic increasing, as it represents time.'
        bch = (baseline_hazards_
                .cumsum()
                .rename('baseline_cumulative_hazards'))
        if set_hazards:
            self.baseline_hazards_ = baseline_hazards_
            self.baseline_cumulative_hazards_ = bch
        return bch

    def predict_cumulative_hazards(self, input, max_duration=None, batch_size=8224, verbose=False,
                                   baseline_hazards_=None, eval_=True, num_workers=0):
        """See `predict_survival_function`."""
        if type(input) is pd.DataFrame:
            input = self.df_to_input(input)
        if baseline_hazards_ is None:
            if not hasattr(self, 'baseline_hazards_'):
                raise ValueError('Need to compute baseline_hazards_. E.g run `model.compute_baseline_hazards()`')
            baseline_hazards_ = self.baseline_hazards_
        assert baseline_hazards_.index.is_monotonic_increasing,\
            'Need index of baseline_hazards_ to be monotonic increasing, as it represents time.'
        return self._predict_cumulative_hazards(input, max_duration, batch_size, verbose, baseline_hazards_,
                                                eval_, num_workers=num_workers)

    def _predict_cumulative_hazards(self, input, max_duration, batch_size, verbose, baseline_hazards_,
                                    eval_=True, num_workers=0):
        raise NotImplementedError

    def predict_surv_df(self, input, max_duration=None, batch_size=8224, verbose=False, baseline_hazards_=None,
                        eval_=True, num_workers=0):
        """Predict survival function for `input`. S(x, t) = exp(-H(x, t))
        Require computed baseline hazards.

        Arguments:
            input {np.array, tensor or tuple} -- Input x passed to net.

        Keyword Arguments:
            max_duration {float} -- Don't compute estimates for duration higher (default: {None})
            batch_size {int} -- Batch size (default: {8224})
            baseline_hazards_ {pd.Series} -- Baseline hazards. If `None` used `model.baseline_hazards_` (default: {None})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            pd.DataFrame -- Survival estimates. One columns for each individual.
        """
        return np.exp(-self.predict_cumulative_hazards(input, max_duration, batch_size, verbose, baseline_hazards_,
                                                       eval_, num_workers))

    def predict_surv(self, input, max_duration=None, batch_size=8224, numpy=None, verbose=False,
                     baseline_hazards_=None, eval_=True, num_workers=0):
        """Predict survival function for `input`. S(x, t) = exp(-H(x, t))
        Require compueted baseline hazards.

        Arguments:
            input {np.array, tensor or tuple} -- Input x passed to net.

        Keyword Arguments:
            max_duration {float} -- Don't compute estimates for duration higher (default: {None})
            batch_size {int} -- Batch size (default: {8224})
            numpy {bool} -- 'False' gives tensor, 'True' gives numpy, and None give same as input
                (default: {None})
            baseline_hazards_ {pd.Series} -- Baseline hazards. If `None` used `model.baseline_hazards_` (default: {None})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            pd.DataFrame -- Survival estimates. One columns for each individual.
        """
        surv = self.predict_surv_df(input, max_duration, batch_size, verbose, baseline_hazards_,
                                    eval_, num_workers)
        surv = torch.from_numpy(surv.values.transpose())
        return tt.utils.array_or_tensor(surv, numpy, input)

    def save_net(self, path, **kwargs):
        """Save self.net and baseline hazards to file.

        Arguments:
            path {str} -- Path to file.
            **kwargs are passed to torch.save

        Returns:
            None
        """
        path, extension = os.path.splitext(path)
        if extension == "":
            extension = '.pt'
        super().save_net(path+extension, **kwargs)
        if hasattr(self, 'baseline_hazards_'):
            self.baseline_hazards_.to_pickle(path+'_blh.pickle')

    def load_net(self, path, **kwargs):
        """Load net and hazards from file.

        Arguments:
            path {str} -- Path to file.
            **kwargs are passed to torch.load

        Returns:
            None
        """
        path, extension = os.path.splitext(path)
        if extension == "":
            extension = '.pt'
        super().load_net(path+extension, **kwargs)
        blh_path = path+'_blh.pickle'
        if os.path.isfile(blh_path):
            self.baseline_hazards_ = pd.read_pickle(blh_path)
            self.baseline_cumulative_hazards_ = self.baseline_hazards_.cumsum()

    def df_to_input(self, df):
        input = df[self.input_cols].values
        return input

class _CoxPHBase(_CoxBase):
    def _compute_baseline_hazards(self, input, df_target, max_duration, batch_size, eval_=True, num_workers=0):
        if max_duration is None:
            max_duration = np.inf
        
        labels = df_target[self.duration_col]

        event_set = np.zeros((labels.unique().size,labels.size)) #Create temporaty set on which we store the event for each development time in batch
        risk_set = np.zeros((labels.unique().size, labels.size)) #Create temporary set on which we store the risk set for each relevant time. 
        unique = labels.unique().tolist()
        unique.sort()

        for i in(range(len(unique))): #create correct exposure and risk set
            event_set[i,:] = (labels == unique[i])
            risk_set[i,:] = (df_target[self.duration_col]>= unique[i]) * (df_target[self.truncation_col] < unique[i])

        pred = np.exp(self.predict(input, batch_size, True, eval_, num_workers=num_workers))

        risk_set = risk_set.dot(pred) - 1/2*(event_set.dot(pred).sum(-1).reshape(-1,1))


        # Here we are computing when expg when there are no events.
        #   Could be made faster, by only computing when there are events.
        return (df_target
                .sort_values(self.duration_col)
                .groupby(self.duration_col)
                .agg({self.event_col: 'sum'})                
                .assign(risk = risk_set)
                .pipe(lambda x: x[self.event_col]/x['risk'])
                .fillna(0.)
                .sort_index(ascending=False)
                .iloc[::-1]
                .loc[lambda x: x.index <= max_duration]
                .rename('baseline_hazards'))

    def _predict_cumulative_hazards(self, input, max_duration, batch_size, verbose, baseline_hazards_,
                                    eval_=True, num_workers=0):
        max_duration = np.inf if max_duration is None else max_duration
        if baseline_hazards_ is self.baseline_hazards_:
            bch = self.baseline_cumulative_hazards_
        else:
            bch = self.compute_baseline_cumulative_hazards(set_hazards=False, 
                                                           baseline_hazards_=baseline_hazards_)
        bch = bch.loc[lambda x: x.index <= max_duration]
        expg = np.exp(self.predict(input, batch_size, True, eval_, num_workers=num_workers)).reshape(1, -1)
        return pd.DataFrame(bch.values.reshape(-1, 1).dot(expg), 
                            index=bch.index)


def cox_ph_loss(log_h: Tensor, durations: Tensor, events: Tensor, truncation: Tensor, tie: str) -> Tensor:
    """Loss for CoxPH model. If data is sorted by descending duration, see `cox_ph_loss_sorted`.

    Dependent on the tie handlng, we create different apporaches for the risk set.
    Be aware that for larger batch sizes, the continuous case can be quite slow.

    At first we create matrixes holding the necesarry values. This approach is taken to limit the calculations being done in for loops. 

    """
    idx = durations.sort(descending=False)[1]
    events = events[idx]
    durations = durations[idx]
    truncation = truncation[idx]
    log_h = log_h[idx].view(-1).mul(events) #"Handles Censoring"

    labels = durations.sort(descending=False)[0].to(torch.double)

    event_set = torch.zeros(labels.unique().shape[0],len(labels.tolist()))
    risk_set = torch.zeros(labels.unique().shape[0],len(labels.tolist()))
    
    unique = labels.to(torch.double).unique().tolist()

    if(tie == 'Continuous'):
        log_like = 0
        for i in(range(len(unique))): 
            event_set =  (labels == unique[i]) * (events==1) #Handles censoring
            risk_set = (log_h.exp().mul((durations>= unique[i]) * (truncation < unique[i]) * (events==1)) ).sum().log() #Handles censoring   
            
            h_events = log_h.mul( (labels == unique[i]) ).sum() 
           # events_summed = event_set.matmul(events) 
            log_like = log_like -h_events.sub(risk_set).sum()
        
        #log_h_tmp = log_h.clone()

        #def risk_set_tmp(time):
        #    return( (log_h_tmp.exp().mul((durations>=time) * (truncation<time)) ).sum() )

        event_set = events

    if(tie == 'Breslow'):
        for i in(range(len(unique))): #This handles ties the breslow way (see https://people.math.aau.dk/~rw/Undervisning/DurationAnalysis/Slides/lektion3.pdf slide 24)
            event_set[i,:] = (labels == unique[i]) * (events==1) #Handles censoring
            risk_set[i,:] = (durations>= unique[i]) * (truncation < unique[i]) * (events==1) #Handles censoring

    
        h_events = event_set.matmul(log_h) #the breslow tie approach (simple sum all covaries that have equal ties)
        h_risk = (risk_set.matmul(log_h.exp()).sub((event_set.matmul(log_h.exp())).mul(1/2))).log() #the 1/2 approach since not all claims should count

        #if no exposure in relevant risk group, set to 0 (would also be handled by events_summed multiplication, but inf*0 = inf in python)
        h_risk[h_risk.abs()==float('inf')]= 0
        h_risk = h_risk.nan_to_num() 

        events_summed = event_set.matmul(events) 
        log_like = -h_events.sub(h_risk.mul(events_summed)).sum()
    
    if(tie == 'Efron'):
        efron_set  = torch.zeros(labels.unique().shape[0],len(labels.tolist()))
        for i in(range(len(unique))): #This handles ties the efron way (see https://people.math.aau.dk/~rw/Undervisning/DurationAnalysis/Slides/lektion3.pdf slide 24)
            event_set[i,:] = (labels == unique[i]) * (events==1) #Handles censoring
            risk_set[i,:] = (durations>= unique[i]) * (truncation < unique[i])
            tmp =  torch.arange(0,torch.sum(event_set[i,:]))
            efron_set[i, tmp.int().tolist()] = tmp #Diffrent set,since we need a special risk set for each claim during a tied event time.


        h_risk = risk_set.matmul(log_h.exp())

        h_events0 = event_set.matmul(log_h.exp())

        efron_risk = (1/torch.sum(event_set,1)).mul(h_events0).reshape(-1,1).matmul(torch.ones(1,len(labels.tolist())))

        efron_risk = efron_risk.mul(efron_set)

        #if no exposure in relevant risk group, set to 0 (would also be handled by events_summed multiplication, but inf*0 = inf in python)
        h_risk[h_risk.abs()==float('inf')]= 1

        events_sort,_ = event_set.sort(-1, descending=True)
        efron_risk = (h_risk.reshape(-1,1).sub(efron_risk)).log().matmul(events_sort.t()).diagonal() #only need the diagonal column,can probably be done smarter
        efron_risk[efron_risk.abs()==float('inf')]= 0

        h_events = event_set.matmul(log_h)


        log_like = -h_events.sub (efron_risk ).sum()
            
    return log_like.div(event_set.sum()) #"return average log likelihood"

class CoxPHLoss(torch.nn.Module):
    """
    Elastic loss defined during init procedure.
    

    """
    def __init__(self, xi=None, eps=None, net=None, tie="Breslow"):
        super().__init__()
        self.xi = xi #for elastic regularization
        self.eps = eps #for elastic regularizaition
        self.net = net
        self.tie = tie #for tie-handling, default = "Breslow"
        
    
    def forward(self, log_h: Tensor, durations: Tensor, events: Tensor, truncation: Tensor) -> Tensor:

        params = []
        for layers, _ in self.net.state_dict().items():
            if layers.find('weight')!=-1:
                 params.append(self.net.get_parameter(layers).view(-1))    
        params = torch.cat(params)
        

        return cox_ph_loss(log_h, durations, events, truncation, self.tie) + self.eps * ( self.xi* (params.square().sum()) + (1-self.xi)*params.abs().sum()) #elastic regularization


class CoxPH(_CoxPHBase):
    """Cox proportional hazards model parameterized with a neural net.
    This is essentially the DeepSurv method [1].

    The loss function is not quite the partial log-likelihood, but close.    
    The difference is that for tied events, we use a random order instead of 
    including all individuals that had an event at that point in time.

    Arguments:
        net {torch.nn.Module} -- A pytorch net.
    
    Keyword Arguments:
        optimizer {torch or torchtuples optimizer} -- Optimizer (default: {None})
        device {str, int, torch.device} -- Device to compute on. (default: {None})
            Preferably pass a torch.device object.
            If 'None': use default gpu if available, else use cpu.
            If 'int': used that gpu: torch.device('cuda:<device>').
            If 'string': string is passed to torch.device('string').

    [1] Jared L. Katzman, Uri Shaham, Alexander Cloninger, Jonathan Bates, Tingting Jiang, and Yuval Kluger.
        Deepsurv: personalized treatment recommender system using a Cox proportional hazards deep neural network.
        BMC Medical Research Methodology, 18(1), 2018.
        https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1
    """
    def __init__(self, net, optimizer=None, device=None, loss=None, xi=None, eps=None, tie="Breslow"):
        if loss is None:
            loss = CoxPHLoss(xi, eps, net, tie)
        super().__init__(net,loss,optimizer,device)