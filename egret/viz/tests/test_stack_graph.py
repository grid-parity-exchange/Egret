import os
import json
import unittest

try:
    import matplotlib
    import seaborn
except (ImportError, ModuleNotFoundError):
    viz_packages_installed = False
else:
    viz_packages_installed = True

    from egret.viz.generate_graphs import generate_stack_graph
    from egret.data.model_data import ModelData
    from egret.models.unit_commitment import solve_unit_commitment, create_tight_unit_commitment_model


@unittest.skipUnless(viz_packages_installed, "matplotlib and seaborn packages are both required to run and test the visualization capabilities.")
class TestStackGraphWithUCTestInstances(unittest.TestCase):
    """Test class for running the unit commitment on the UC test instances and generating stack graphs for each."""
    def setUp(self):
        pass
    
    def test_standard_stack_graph(self):
        """Tests standard stack graph generation."""
        current_dir = os.path.dirname(os.path.abspath(__file__))
        test_cases = [os.path.join(current_dir, '..', '..', 'models', 'tests', 'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1,6)]

        for test_case in test_cases:   
            with open(test_case, 'r') as f:
                md_dict = json.load(f)
            md = ModelData(md_dict)

            solved_md = solve_unit_commitment(md,
                            'cbc',
                            mipgap = 0.001,
                            timelimit = None,
                            solver_tee = True,
                            symbolic_solver_labels = False,
                            solver_options = None,
                            uc_model_generator=create_tight_unit_commitment_model,
                            relaxed=False,
                            return_model=False)

            fig, ax = generate_stack_graph(
                solved_md, 
                title=repr(test_case),
                show_individual_components=False,
                plot_individual_generators=False,
            )
    
    def test_individual_component_stack_graph(self):
        """Tests stack graph generation when breaking out individual components per generation type."""
        current_dir = os.path.dirname(os.path.abspath(__file__))
        test_cases = [os.path.join(current_dir, '..', '..', 'models', 'tests', 'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1, 2)]

        for test_case in test_cases:   
            with open(test_case, 'r') as f:
                md_dict = json.load(f)
            md = ModelData(md_dict)

            solved_md = solve_unit_commitment(md,
                            'cbc',
                            mipgap = 0.001,
                            timelimit = None,
                            solver_tee = True,
                            symbolic_solver_labels = False,
                            solver_options = None,
                            uc_model_generator=create_tight_unit_commitment_model,
                            relaxed=False,
                            return_model=False)

            fig, ax = generate_stack_graph(
                solved_md, 
                title=repr(test_case),
                show_individual_components=True,
                plot_individual_generators=False,
            )
    
    def test_individual_generator_stack_graph(self):
        """Tests stack graph generation when plotting individual generators."""
        current_dir = os.path.dirname(os.path.abspath(__file__))
        test_cases = [os.path.join(current_dir, '..', '..', 'models', 'tests', 'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1, 2)]

        for test_case in test_cases:   
            with open(test_case, 'r') as f:
                md_dict = json.load(f)
            md = ModelData(md_dict)

            solved_md = solve_unit_commitment(md,
                            'cbc',
                            mipgap = 0.001,
                            timelimit = None,
                            solver_tee = True,
                            symbolic_solver_labels = False,
                            solver_options = None,
                            uc_model_generator=create_tight_unit_commitment_model,
                            relaxed=False,
                            return_model=False)

            with self.assertRaises(ValueError):
                # The number of generators in the test case exceeds the maximum for this feature.
                fig, ax = generate_stack_graph(
                    solved_md, 
                    title=repr(test_case),
                    show_individual_components=False,
                    plot_individual_generators=True,
                )
    
    def test_individual_generator_stack_graph_exception(self):
        """Tests for stack graph generation to fail when both 'show_individual_components' and 'plot_individual_generators' are simultaneously True."""
        current_dir = os.path.dirname(os.path.abspath(__file__))
        test_cases = [os.path.join(current_dir, '..', '..', 'models', 'tests', 'uc_test_instances', 'test_case_{}.json'.format(i)) for i in range(1, 2)]

        for test_case in test_cases:   
            with open(test_case, 'r') as f:
                md_dict = json.load(f)
            md = ModelData(md_dict)

            solved_md = solve_unit_commitment(md,
                            'cbc',
                            mipgap = 0.001,
                            timelimit = None,
                            solver_tee = True,
                            symbolic_solver_labels = False,
                            solver_options = None,
                            uc_model_generator=create_tight_unit_commitment_model,
                            relaxed=False,
                            return_model=False)

            with self.assertRaises(ValueError):
                # You cannot set both options to True simultaneously.
                fig, ax = generate_stack_graph(
                    solved_md, 
                    title=repr(test_case),
                    show_individual_components=True,
                    plot_individual_generators=True,
                )


if __name__ == '__main__':
    unittest.main()
