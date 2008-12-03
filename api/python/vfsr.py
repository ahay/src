import c_vfsr

class VFSR_defines:
    def __init__(self):
        self.this = c_vfsr.VFSR_DEFINES()
        self.this.COST_PRECISION = 1.0E-18
        self.this.USER_INITIAL_PARAMETERS = 0
        self.this.ACCEPTED_TO_GENERATED_RATIO = 1.0E-4
        self.this.LIMIT_ACCEPTANCES = 1000
        self.this.TEMPERATURE_RATIO_SCALE = 1.0E-5
        self.this.TEMPERATURE_ANNEAL_SCALE = 100.0
        self.this.COST_PARAMETER_SCALE = 1.0
        self.this.TESTING_FREQUENCY_MODULUS = 100
        self.this.MAXIMUM_REANNEAL_INDEX = 100000000
        self.this.REANNEAL_RESCALE = 10.0
        self.this.INITIAL_PARAMETER_TEMPERATURE = 1.0
        self.this.USER_INITIAL_PARAMETERS_TEMPS = 0
        self.this.USER_INITIAL_COST_TEMP = 0
        self.this.NUMBER_COST_SAMPLES = 5
        self.this.MAXIMUM_COST_REPEAT = 5
        self.this.DELTA_X = 0.001
        self.this.INCLUDE_INTEGER_PARAMETERS = 0
        self.this.ACTIVATE_REANNEAL = 1
        self.this.LIMIT_INVALID_GENERATED_STATES = 1000
    def cost_precision(self,val):
        self.this.COST_PRECISION = val
    def user_initial_parameters(self,val):
        self.this.USER_INITIAL_PARAMETERS = val
    def accepted_to_generated_ratio(self,val):
        self.this.ACCEPTED_TO_GENERATED_RATIO = val
    def limit_acceptances(self,val):
        self.this.LIMIT_ACCEPTANCES = 10000
    def temperature_ratio_scale(self,val):
        self.this.TEMPERATURE_RATIO_SCALE = val
    def temperature_anneal_scale(self,val):
        self.this.TEMPERATURE_ANNEAL_SCALE = val
    def cost_parameter_scale(self,val):
        self.this.COST_PARAMETER_SCALE = val
    def testing_frequency_modulus(self,val):
        self.this.TESTING_FREQUENCY_MODULUS = val
    def maximum_reanneal_index(self,val):
        self.this.MAXIMUM_REANNEAL_INDEX = val
    def reanneal_scale(self,val):
        self.this.REANNEAL_RESCALE = val
    def initial_parameter_temperature(self,val):
        self.this.INITIAL_PARAMETER_TEMPERATURE = val
    def user_initial_parameters_temps(self,val):
        self.this.USER_INITIAL_PARAMETERS_TEMPS = val
    def user_initial_cost_temp(self,val):
        self.this.USER_INITIAL_COST_TEMP = val
    def number_cost_samples(self,val):
        self.this.NUMBER_COST_SAMPLES = val
    def maximum_cost_repeat(self,val):
        self.this.MAXIMUM_COST_REPEAT = val
    def delta_x(self,val):
        self.this.DELTA_X = val
    def include_integer_parameters(self,val):
        self.this.INCLUDE_INTEGER_PARAMETERS = 0
    def activate_reanneal(self,val):
        self.this.ACTIVATE_REANNEAL = val
    def limit_invalid_generated_states(self,val):
        self.this.LIMIT_INVALID_GENERATED_STATES = val


def run(cost_func, parameter_type,
        parameter_initial_final,
        parameter_minimum, parameter_maximum,
        tangents, curvature,
        exit_status,
        OPTIONS):
    return c_vfsr.run(cost_func,
                      parameter_type, parameter_initial_final,
                      parameter_minimum, parameter_maximum,
                      tangents, curvature,
                      exit_status, OPTIONS.this)
