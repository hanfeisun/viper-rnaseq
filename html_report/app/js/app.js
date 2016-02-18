'use strict';

/* App Module */

var viperApp = angular.module('viperApp', [
  'ui.router',
  'viperAppControllers'
]);


viperApp.config(['$stateProvider', '$urlRouterProvider', '$locationProvider', 
        function($stateProvider, $urlRouterProvider, $locationProvider) {

    $locationProvider.html5Mode({
        enabled: true,
        requireBase: false
    });    

    $urlRouterProvider.otherwise('/');
    var state_list = ["Alignment-Report", "Read-Distribution-Report", "PCA-Plots", "Heatmaps", "Volcano-Plots"];
    for(var index=0; index < state_list.length; index++) {
        $stateProvider.state('/' + state_list[index], {
            url: '/' + state_list[index],
            templateUrl: '/app/partials/section-report.html'
        });
    }
    $stateProvider

        .state('/', {
            url: '/',
            templateUrl: '/app/partials/full-report.html'
        })
        .state('/home', {
            url: '/home',
            templateUrl: '/app/partials/home.html'
        })
        .state('/about',{
            url: '/about',
            templateUrl: '/app/partials/about.html'
        })
        .state('/contact',{
            url: '/contact',
            templateUrl: '/app/partials/contact.html'
        });
}]);


